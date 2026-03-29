# -*- coding: utf-8 -*-
"""
OSMgen: Nástroj pro hybridní topologický model a generalizaci říční sítě
========================================================================
Nástroj vyvinutý pro extrakci, topologickou analýzu a kartografickou 
generalizaci říčních sítí ze surových dat OpenStreetMap (OSM).
Model je optimalizován pro přehledná mapová měřítka.

Klíčové algoritmy:
- Tuple-Based Scout (Rozhodování na bifurkacích)
- Kahnův topologický algoritmus (Třídění pro akumulaci)
- Zpětná propagace (Ochrana proti pahýlům a ztrátě pramenů)
- BFS Ink Propagation (Záchrana paralelních ramen a ostrovů)
- Pairwise Buffer (Rychlé generování plošného koryta)
"""

import arcpy
import math
import os
import unicodedata
import sys
from collections import defaultdict, deque
from datetime import datetime

# =========================================================
# 1. POMOCNÉ FUNKCE A TEXTOVÁ NORMALIZACE
# =========================================================

def log(msg):
    """Zápis zpráv do konzole ArcGIS a standardního výstupu."""
    arcpy.AddMessage(msg)
    print(msg)

def normalize_name(name):
    """
    Odstraní diakritiku a převede text na malá písmena.
    Zajišťuje bezpečné a jednotné porovnávání textových atributů.
    """
    if not name: return ""
    return "".join(c for c in unicodedata.normalize('NFD', str(name).lower()) if unicodedata.category(c) != 'Mn').strip()

def get_name_words(name):
    """
    Rozloží název na množinu jednotlivých slov pro částečné shody 
    (např. 'Černá Ostravice' -> {'cerna', 'ostravice'}).
    """
    norm = normalize_name(name)
    return set(norm.replace('-', ' ').replace(',', ' ').split()) if norm else set()

def scout_branch(source_oid, target_oid, river_data, downstream_targets, target_name_words, max_depth=8):
    """
    Heuristický algoritmus prohledávání (Lookahead).
    Hodnotí vhodnost cesty na soutocích nahlížením až o `max_depth` segmentů vpřed.
    Vrací multikriteriální skóre (Tuple), kde má prioritu shoda názvu nad délkou.
    """
    paths = []
    stack = [[target_oid]]
    
    # Průchod grafem do definované hloubky
    while stack:
        path = stack.pop()
        curr = path[-1]
        
        if len(path) == max_depth:
            paths.append(path)
            continue
            
        nxt_opts = [t for t in downstream_targets.get(curr, []) if t in river_data and t not in path]
        if not nxt_opts:
            paths.append(path)
        else:
            for n in nxt_opts:
                stack.append(path + [n])
                
    best_score = (-1, -1, -1, -1) 
    
    # Vyhodnocení nalezených tras
    for path in paths:
        exact_names = sum(1 for n in path if target_name_words and river_data[n]["words"] == target_name_words)
        partial_names = sum(1 for n in path if target_name_words and river_data[n]["words"].intersection(target_name_words))
        rivers = sum(1 for n in path if river_data[n]["fclass"] == 'river')
        total_len = sum(river_data[n]["len"] for n in path)
        
        # Tuple skórování: Python automaticky hodnotí prioritu zleva doprava
        score = (exact_names, partial_names, rivers, total_len)
        if score > best_score:
            best_score = score
            
    return best_score

def get_param(index, default_value):
    """Bezpečné načtení parametrů z GUI nástroje ArcGIS."""
    try:
        val = arcpy.GetParameterAsText(index)
        if not val: return default_value
        if isinstance(default_value, float): return float(val.replace(',', '.'))
        if isinstance(default_value, int): return int(float(val.replace(',', '.')))
        return val
    except:
        return default_value

# ==============================================================
# 2. INICIALIZACE PARAMETRŮ A PRACOVNÍHO PROSTŘEDÍ
# ==============================================================

SRC_LINES       = arcpy.GetParameterAsText(0) or r"C:\Users\adams\Documents\BP_project_GIS\DATA IN\gis_osm_waterways_free_1.shp"
TARGET_SCALE    = get_param(1, 100000.0)
MAIN_RIVER_NAME = arcpy.GetParameterAsText(2)
OSM_ID_INPUT    = arcpy.GetParameterAsText(3)
OUT_GDB         = arcpy.GetParameterAsText(4) 
BASE_RADIUS     = get_param(5, 0.32) 

# Načtení cílového souřadnicového systému (jako prostorového objektu)
try:
    OUT_SR = arcpy.GetParameter(6)
except:
    OUT_SR = None

ts = datetime.now().strftime("%Y%m%d_%H%M%S")
raw_name = normalize_name(MAIN_RIVER_NAME)
safe_river = "".join(c if c.isalnum() else "_" for c in raw_name).strip("_").capitalize() if raw_name else "Reka"

# Priorita názvu výstupů podle zadaných parametrů
if OSM_ID_INPUT: safe_river = f"OSM_{str(OSM_ID_INPUT).strip()}"

# Zajištění výstupní geodatabáze
if OUT_GDB and os.path.exists(OUT_GDB):
    gdb_path = OUT_GDB
else:
    current_ws = arcpy.env.workspace
    if current_ws and current_ws.lower().endswith('.gdb'): 
        OUT_FOLDER = os.path.dirname(current_ws)
    else: 
        OUT_FOLDER = current_ws if current_ws else arcpy.env.scratchFolder

    if not os.path.exists(OUT_FOLDER):
        try: os.makedirs(OUT_FOLDER)
        except: pass
    
    gdb_name = f"Povodi_{safe_river}_{ts}"
    gdb_path = os.path.join(OUT_FOLDER, f"{gdb_name}.gdb")
    arcpy.management.CreateFileGDB(OUT_FOLDER, gdb_name)

arcpy.env.workspace = gdb_path

# Definiční knihovny pro sémantickou analýzu
cutoff_raw       = "umělý kanál, přivaděč, náhon, kanál"
CUTOFF_NAMES     = [w for n in cutoff_raw.split(',') for w in get_name_words(n)]
WATERBODY_KW     = {"nadrz", "rybnik", "jezero", "prehrada", "tajch", "ryb.", "vodni"}
SUSPICIOUS_LIMIT = 50000.0 
NAME_FIELD       = "name"

arcpy.env.overwriteOutput = True
arcpy.env.parallelProcessingFactor = "100%"

log(f"=== OSMgen: SPUŠTĚNO PRO CÍLOVÉ MĚŘÍTKO 1:{int(TARGET_SCALE)} ===")
log(f"=== VÝSTUPNÍ GEODATABÁZE: {gdb_path} ===")

# =========================================================
# FÁZE A: GEOMETRICKÁ DEKOMPOZICE A PROJEKCE
# =========================================================

log(f"\n[FÁZE A] 1) Geometrická dekompozice linií na soutocích...")
temp_planar_name = arcpy.ValidateTableName(f"L01_Planarized_{ts}", gdb_path)
temp_planar = os.path.join(gdb_path, temp_planar_name)

# Rozdělení multipart prvků a vynucení topologických uzlů (ještě v původní projekci)
arcpy.management.FeatureToLine(SRC_LINES, temp_planar, attributes="ATTRIBUTES")

log(f"[FÁZE A] 2) Transformace do metrického souřadnicového systému...")
input_desc = arcpy.Describe(temp_planar)
spatial_ref = input_desc.spatialReference

# Ošetření nezadaného systému (fallback na národní systém S-JTSK)
if OUT_SR and getattr(OUT_SR, 'name', None):
    target_sr = OUT_SR
else:
    target_sr = arcpy.SpatialReference(5514) 

lines_fc_name = arcpy.ValidateTableName(f"L02_WorkingLines_{ts}", gdb_path)
lines_fc = os.path.join(gdb_path, lines_fc_name)

if spatial_ref.factoryCode != target_sr.factoryCode:
    log(f"  -> Aplikuji transformaci z {spatial_ref.name} do {target_sr.name}...")
    arcpy.management.Project(temp_planar, lines_fc, target_sr)
else:
    arcpy.management.CopyFeatures(temp_planar, lines_fc)

arcpy.management.Delete(temp_planar)

# =========================================================
# FÁZE B: EXTRAKCE TOPOLOGIE A TRASOVÁNÍ GRAFU
# ========================================================= 

log(f"\n[FÁZE B] 1) Budování struktur orientovaného grafu (DAG)...")
end_pts_name = arcpy.ValidateTableName(f"temp_end_{ts}", gdb_path)
start_pts_name = arcpy.ValidateTableName(f"temp_start_{ts}", gdb_path)
temp_end = os.path.join(gdb_path, end_pts_name)
temp_start = os.path.join(gdb_path, start_pts_name)

# Extrakce koncových bodů pro definici směru proudění
arcpy.management.FeatureVerticesToPoints(lines_fc, temp_end, "END")
arcpy.management.FeatureVerticesToPoints(lines_fc, temp_start, "START")

near_name = arcpy.ValidateTableName(f"temp_near_{ts}", gdb_path)
temp_near = os.path.join(gdb_path, near_name)

# Topologický magnet: Zvýšená tolerance (1.5 m) přemosťuje drobné datové chyby (mikromezery)
arcpy.analysis.GenerateNearTable(temp_end, temp_start, temp_near, "1.5 Meters", closest="ALL")

endpt_to_seg = {r[0]: r[1] for r in arcpy.da.SearchCursor(temp_end, ["OID@", "ORIG_FID"])}
startpt_to_seg = {r[0]: r[1] for r in arcpy.da.SearchCursor(temp_start, ["OID@", "ORIG_FID"])}

# Načtení atributů sítě do paměti
river_data = {}
fields_in_fc = [f.name for f in arcpy.ListFields(lines_fc)]
valid_class_fields = ['fclass', 'waterway', 'type', 'class']
class_field = next((f for f in fields_in_fc if f.lower() in valid_class_fields), None)
osm_id_field = next((f for f in fields_in_fc if f.lower() == 'osm_id'), None)

search_fields = ["OID@", NAME_FIELD, "SHAPE@LENGTH", "SHAPE@"]
if class_field: search_fields.append(class_field)
if osm_id_field: search_fields.append(osm_id_field)

with arcpy.da.SearchCursor(lines_fc, search_fields) as cur:
    for row in cur:
        oid, nm, length, geom = row[0], row[1], row[2], row[3]
        fcls = str(row[search_fields.index(class_field)]).lower() if class_field else ""
        row_osm_id = str(row[search_fields.index(osm_id_field)]) if osm_id_field and row[search_fields.index(osm_id_field)] else ""
        
        river_data[oid] = {
            "name": nm if nm else "", "norm": normalize_name(nm), "words": get_name_words(nm), 
            "fclass": fcls, "len": length, "osm_id": row_osm_id
        }

downstream_targets = defaultdict(list)
upstream_sources = defaultdict(list)

# Sestavení paměťových seznamů sousednosti (Adjacency Lists)
with arcpy.da.SearchCursor(temp_near, ["IN_FID", "NEAR_FID"]) as cur:
    for pt_end, pt_start in cur:
        source_seg = endpt_to_seg.get(pt_end)
        target_seg = startpt_to_seg.get(pt_start)
        if source_seg and target_seg and source_seg != target_seg:
            if source_seg in river_data and target_seg in river_data:
                downstream_targets[source_seg].append(target_seg)
                upstream_sources[target_seg].append(source_seg)

log("\n[FÁZE B] 2) Definice zájmového toku a lokalizace kořenového uzlu...")
start_node = None
target_name_words = set()

# Primární pokus: Vyhledání podle OSM ID
if OSM_ID_INPUT:
    candidates = [oid for oid, d in river_data.items() if d["osm_id"] == str(OSM_ID_INPUT).strip()]
    if candidates:
        start_node = max(candidates, key=lambda x: river_data[x]["len"])
        target_name_words = river_data[start_node]["words"]
        log(f"  -> Uzel inicializován na základě OSM ID ({OSM_ID_INPUT}).")

# Sekundární pokus: Vyhledání podle textového názvu toku
if not start_node and MAIN_RIVER_NAME:
    main_river_words = get_name_words(MAIN_RIVER_NAME.replace(',', ' '))
    candidates = [oid for oid, d in river_data.items() if d["words"].intersection(main_river_words)]
    if candidates:
        start_node = max(candidates, key=lambda x: river_data[x]["len"])
        target_name_words = main_river_words
        log(f"  -> Uzel inicializován na základě sémantické shody.")

if not start_node:
    arcpy.AddError("Kritická chyba: Počáteční segment nebyl v datech nalezen.")
    sys.exit()

log("\n[FÁZE B] 3) Evaluace bifurkací (Scout Algoritmus)...")
primary_downstream = {}
primary_upstream = defaultdict(list)

# Stanovení primární větve na každém větvení
for source, targets in downstream_targets.items():
    valid_targets = [t for t in targets if t in river_data]
    if not valid_targets: continue
    
    best_target = max(valid_targets, key=lambda t: scout_branch(source, t, river_data, downstream_targets, target_name_words, max_depth=8))
    primary_downstream[source] = best_target
    primary_upstream[best_target].append(source)

log("\n[FÁZE B] 4) Detekce ústí a vymezení hydrologického povodí (DFS)...")
curr = start_node
visited_down = set([curr])

# Walk-down fáze: Trasování k ústí řeky
while True:
    nxt = primary_downstream.get(curr)
    if not nxt or nxt in visited_down: break
    curr_words = river_data[curr]["words"]
    nxt_words = river_data[nxt]["words"]
    other_tributaries = [t for t in upstream_sources[nxt] if t != curr]
    is_demise = False
    
    # Detekce zániku toku s ohledem na říční ostrovy a paralelní ramena
    if nxt_words:
        we_continue = bool(curr_words.intersection(nxt_words))
        for t in other_tributaries:
            if bool(river_data[t]["words"].intersection(nxt_words)) and not we_continue:
                is_demise = True; break
                
    if not is_demise and target_name_words and nxt_words:
        if not target_name_words.intersection(nxt_words) and not WATERBODY_KW.intersection(nxt_words): is_demise = True

    if is_demise: break
    visited_down.add(nxt)
    curr = nxt
    
outlet_oid = curr

# DFS fáze: Sběr všech nadřazených větví (izolace povodí)
stack = [outlet_oid]
basin_nodes = set([outlet_oid])
while stack:
    curr = stack.pop()
    for c in primary_upstream[curr]:
        if c not in basin_nodes:
            basin_nodes.add(c); stack.append(c)

# =========================================================
# FÁZE C: VÝPOČET HYDROMETRIK A OCHRANNÝCH FILTRŮ
# =========================================================

log("\n[FÁZE C] 1) Topologické třídění sítě (Kahnův algoritmus)...")
for field in ["upstr_len", "shreve_proxy", "max_ds_shreve", "filter_val"]: arcpy.management.AddField(lines_fc, field, "DOUBLE")
for field in ["shreve_val", "strahler_val"]: arcpy.management.AddField(lines_fc, field, "LONG")
arcpy.management.AddField(lines_fc, "is_basin", "SHORT")

in_degree = {n: 0 for n in basin_nodes}
for n in basin_nodes:
    for c in primary_upstream[n]:
        if c in basin_nodes: in_degree[n] += 1

# Řazení uzlů od pramenů k ústí
queue = deque([n for n in basin_nodes if in_degree[n] == 0])
topo_order, remaining = [], set(basin_nodes)

while remaining:
    if not queue:
        forced_node = next(iter(remaining))
        in_degree[forced_node] = 0; queue.append(forced_node)
        
    curr = queue.popleft()
    if curr not in remaining: continue
    remaining.remove(curr); topo_order.append(curr)
        
    target = primary_downstream.get(curr)
    if target and target in remaining:
        in_degree[target] -= 1
        if in_degree[target] == 0: queue.append(target)

log("[FÁZE C] 2) Výpočet Strahlerova a Shreveho řádu toku...")
upstr_len, shreve_vals, strahler_vals = {}, {}, {}

for node in topo_order:
    primary_children = [c for c in primary_upstream[node] if c in basin_nodes and c in upstr_len]
    if not primary_children:
        upstr_len[node], shreve_vals[node], strahler_vals[node] = river_data[node]["len"], 1, 1
    else:
        upstr_len[node] = river_data[node]["len"] + sum(upstr_len[c] for c in primary_children)
        shreve_vals[node] = sum(shreve_vals[c] for c in primary_children)
        c_str = [strahler_vals[c] for c in primary_children]
        max_s = max(c_str)
        strahler_vals[node] = max_s + 1 if c_str.count(max_s) >= 2 else max_s

log("[FÁZE C] 3) Ochrana pramenů: Zpětná propagace (Inherited Shreve a Length)...")
max_ds_shreve = {}
max_ds_len = {} # Ochrana před pahýly
main_trunk_nodes = set([outlet_oid]) 

for node in reversed(topo_order):
    if node not in max_ds_shreve:
        max_ds_shreve[node] = shreve_vals[node]
    if node not in max_ds_len:
        max_ds_len[node] = upstr_len[node]
        
    primary_children = [c for c in primary_upstream[node] if c in basin_nodes and c in upstr_len]
    if primary_children:
        # Volba hlavní větve pro propagaci štítu
        main_child = max(primary_children, key=lambda c: (
            bool(target_name_words and river_data[c]["words"].intersection(target_name_words)),
            shreve_vals[c], 
            upstr_len[c]
        ))
        
        for c in primary_children:
            # Hlavní dítě dědí hodnoty zespodu, vedlejší děti začínají vlastní akumulaci
            max_ds_shreve[c] = max_ds_shreve[node] if c == main_child else shreve_vals[c]
            max_ds_len[c] = max_ds_len[node] if c == main_child else upstr_len[c]
            
            if node in main_trunk_nodes and c == main_child:
                main_trunk_nodes.add(c)

log("[FÁZE C] 4) Inicializace sémantického štítu a délkových filtrů...")
oid_proxy, oid_filters = {}, {}

# Vyloučení generických slov (aby každý "potok" nedostal imunitu proti smazání)
generic_words = {"potok", "reka", "r.", "p.", "kanal", "nadrz", "vodni"}
specific_name_words = target_name_words - generic_words
if not specific_name_words: specific_name_words = target_name_words

# Exponenciální prahy pro čistší malá měřítka
LEN_THRESHOLD = TARGET_SCALE * 0.12  
SHR_THRESHOLD = max(1.0, (TARGET_SCALE / 20000.0) ** 2) 

for node in topo_order:
    l_val = upstr_len[node]
    sh_val = shreve_vals[node]
    inherited_shreve = max_ds_shreve.get(node, sh_val)
    inherited_len = max_ds_len.get(node, l_val) # TADY JE TA MAGIE PROTI PAHÝLŮM
    is_true_main = (node in main_trunk_nodes)
    
    shreve_proxy = math.sqrt(max(1.0, sh_val))
    
    if is_true_main:
        f_val = 1.0 
    else:
        # Hodnotíme CELÝ tok pomocí inherited_len, neřežeme ho na kousky!
        f_val = 1.0 if (inherited_len >= LEN_THRESHOLD or inherited_shreve >= SHR_THRESHOLD) else 0.0
        
        # Eliminace umělých kanálů
        if any(c in river_data[node]["words"] for c in CUTOFF_NAMES) and l_val > SUSPICIOUS_LIMIT:
            shreve_proxy, f_val = 0.0, 0.0

    oid_proxy[node], oid_filters[node] = shreve_proxy, f_val

# Zajištění imunity pro Topologickou páteř (na základě specifických jmen)
backbone = set()
if target_name_words:
    curr = outlet_oid
    backbone.add(curr)
    while curr:
        parents = [p for p in primary_upstream.get(curr, []) if p in basin_nodes]
        if not parents: break
        best_p = max(parents, key=lambda p: (bool(river_data[p]["words"].intersection(specific_name_words)), river_data[p]["len"]))
        backbone.add(best_p)
        curr = best_p
        
    for node in basin_nodes:
        is_same_river = bool(river_data[node]["words"].intersection(specific_name_words))
        if node in backbone or is_same_river:
            oid_filters[node] = 1.0

log("[FÁZE C] 5) Ochrana ostrovů: Prohledávání do šířky (Ink Propagation)...")
ink_queue = deque()
for node in topo_order:
    alfa_tgt = primary_downstream.get(node)
    for t in [tgt for tgt in downstream_targets.get(node, []) if tgt in basin_nodes]:
        if t != alfa_tgt:
            # Propustí část vizuální váhy do vedlejších ramen ostrovů
            is_named_branch = bool(specific_name_words and river_data[t]["words"].intersection(specific_name_words))
            multiplier = 0.85 if is_named_branch else 0.15
            ink_queue.append((t, oid_proxy.get(node, 0) * multiplier, oid_filters.get(node, 0)))

while ink_queue:
    curr, proxy, flt = ink_queue.popleft()
    
    # Ochrana proti přepsání hlavní páteře špatnou šířkou z ostrova
    is_true_main = (curr in main_trunk_nodes)
    
    if not is_true_main:
        updated = False
        
        # 1. Nezávislé šíření imunity (Filtru) - ZACHRÁNÍ OSTROV PŘED SMAZÁNÍM
        if flt > oid_filters.get(curr, 0):
            oid_filters[curr] = flt
            updated = True
            
        # 2. Nezávislé šíření vizuální šířky (Proxy) - UDĚLÁ OSTROV UŽŠÍ NEŽ HLAVNÍ TOK
        if proxy > oid_proxy.get(curr, -1):
            oid_proxy[curr] = proxy
            updated = True
            
        if updated:
            for t in [tgt for tgt in downstream_targets.get(curr, []) if tgt in basin_nodes]:
                ink_queue.append((t, proxy, flt)) 

log("[FÁZE C] 6) Zápis vypočtených metrik do atributové tabulky (DEBUG vrstva)...")
with arcpy.da.UpdateCursor(lines_fc, ["OID@", "upstr_len", "shreve_val", "strahler_val", "shreve_proxy", "max_ds_shreve", "filter_val", "is_basin"]) as cur:
    for row in cur:
        oid = row[0]
        if oid in basin_nodes:
            row[1], row[2], row[3] = upstr_len.get(oid, 0), shreve_vals.get(oid, 0), strahler_vals.get(oid, 0)
            row[4], row[5], row[6] = oid_proxy.get(oid, 0), max_ds_shreve.get(oid, 0), oid_filters.get(oid, 0)
            row[7] = 1 
            cur.updateRow(row)
        else:
            row[1], row[2], row[3], row[4], row[5], row[6], row[7] = 0, 0, 0, 0.0, 0, 0.0, 0
            cur.updateRow(row)

debug_name = arcpy.ValidateTableName(f"DEBUG_Sit_{safe_river}_{ts}", gdb_path)
debug_fc = os.path.join(gdb_path, debug_name)
arcpy.management.CopyFeatures(lines_fc, debug_fc)

# =========================================================
# FÁZE D: KARTOGRAFICKÁ GENERALIZACE A 2D KORYTO
# =========================================================

log(f"\n[FÁZE D] Datová redukce a generování polygonů (Cílové měřítko 1:{int(TARGET_SCALE)})...")

sel_name = arcpy.ValidateTableName(f"Temp_Selected_{ts}", gdb_path)
temp_selected_fc = os.path.join(gdb_path, sel_name)
# Odstranění segmentů, které nesplnily podmínky ochrany
arcpy.analysis.Select(lines_fc, temp_selected_fc, "filter_val >= 0.5")

dp_tol = max(10, int(TARGET_SCALE / 2000)) 
log(f"  -> Krok 1: Redukce lomových bodů (Douglas-Peucker alg., tolerance: {dp_tol} m)")
simp_name = arcpy.ValidateTableName(f"Temp_Simplified_{ts}", gdb_path)
temp_simplified_fc = os.path.join(gdb_path, simp_name)
arcpy.cartography.SimplifyLine(temp_selected_fc, temp_simplified_fc, "BEND_SIMPLIFY", f"{dp_tol} Meters", "RESOLVE_ERRORS", "KEEP_COLLAPSED_POINTS")

smooth_tol = max(50, int(TARGET_SCALE / 500)) 
log(f"  -> Krok 2: Vyhlazení geometrie meandrů (PAEK alg., tolerance: {smooth_tol} m)")
osa_name = arcpy.ValidateTableName(f"Osa_{safe_river}_{int(TARGET_SCALE)}_{ts}", gdb_path)
final_carto_fc = os.path.join(gdb_path, osa_name)

# Exaktní zachování pozice krajních bodů (FIXED_CLOSED_ENDPOINT) brání rozpadu sítě
arcpy.cartography.SmoothLine(
    in_features=temp_simplified_fc, 
    out_feature_class=final_carto_fc, 
    algorithm="PAEK", 
    tolerance=f"{smooth_tol} Meters", 
    endpoint_option="FIXED_CLOSED_ENDPOINT", 
    error_option="RESOLVE_ERRORS"
)

log("  -> Inicializace matematické definice poloměru koryta (W = k * sqrt(Shreve))...")
arcpy.management.AddField(final_carto_fc, "buff_dist", "DOUBLE")
with arcpy.da.UpdateCursor(final_carto_fc, ["shreve_proxy", "buff_dist"]) as cur:
    for row in cur:
        row[1] = BASE_RADIUS * row[0] 
        cur.updateRow(row)

log("  -> Tvorba celistvého 2D polygonu (Turbo Pairwise Buffer)...")
poly_name = arcpy.ValidateTableName(f"Koryto_Polygon_{safe_river}_{int(TARGET_SCALE)}_{ts}", gdb_path)
buffer_fc = os.path.join(gdb_path, poly_name)

# Paralelizovaný algoritmus pro masivní úsporu výpočetního času
arcpy.analysis.PairwiseBuffer(final_carto_fc, buffer_fc, "buff_dist", "ALL")

# Průběžný úklid dočasných vrstev
try: 
    for tmp in [temp_end, temp_start, temp_near, lines_fc, temp_selected_fc, temp_simplified_fc]: 
        arcpy.management.Delete(tmp)
except: pass

log("\n" + "="*70)
log(f"  VÝPOČET ÚSPĚŠNĚ DOKONČEN (Souřadnicový systém: {target_sr.name})")
log(f"  2D Koryto uloženo: {buffer_fc}")
log(f"  1D Osa uložena: {final_carto_fc}")
log("="*70)

# Automatické přidání vrstev do aktivní mapy v prostředí ArcGIS Pro
try:
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    m = aprx.activeMap
    if m: 
        m.addDataFromPath(buffer_fc)
        m.addDataFromPath(final_carto_fc)
        m.addDataFromPath(debug_fc)
except Exception: pass
