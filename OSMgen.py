# -*- coding: utf-8 -*-
"""
BAKALÁŘKA 2026: HYBRIDNÍ TOPOLOGICKÝ MODEL ŘÍČNÍ SÍTĚ (BULLETPROOF)
=====================================================
- STRIKTNĚ ORIENTOVANÝ GRAF: Zabraňuje toku proti proudu.
- SCOUT ALGORITMUS: Hledá hlavní větve (Jméno -> Fclass -> Délka).
- RYCHLÝ KAHN: O(1) Cycle Breaker.
- INK PROPAGATION 85/15: Vizuální barvení ostrovů (Nerozbíjí topologii!).
- KARTOGRAFICKÝ INDEX: sqrt(Shreve * Upstream_Length).
- PSEUDO-PROPORTIONAL SYMBOLOGY: 100 tříd (0.1 - 7.0 pt).
"""

import arcpy
import math
import os
import unicodedata
from collections import defaultdict, deque
from datetime import datetime

# =========================================================
# POMOCNÉ FUNKCE
# =========================================================
def log(msg):
    arcpy.AddMessage(msg)
    print(msg)

def normalize_name(name):
    if not name: return ""
    return "".join(c for c in unicodedata.normalize('NFD', str(name).lower()) if unicodedata.category(c) != 'Mn').strip()

def get_name_words(name):
    norm = normalize_name(name)
    return set(norm.replace('-', ' ').replace(',', ' ').split()) if norm else set()

def calc_deviation(naj, geom_out):
    if not geom_in or not geom_out: return 90.0 
    
    dx1 = geom_in['last'][0] - geom_in['first'][0]
    dy1 = geom_in['last'][1] - geom_in['first'][1]
    dx2 = geom_out['last'][0] - geom_out['first'][0]
    dy2 = geom_out['last'][1] - geom_out['first'][1]
    
    mag1 = math.hypot(dx1, dy1)
    mag2 = math.hypot(dx2, dy2)
    
    if mag1 == 0 or mag2 == 0: return 90.0 
    
    cos_val = (dx1 * dx2 + dy1 * dy2) / (mag1 * mag2)
    cos_val = max(-1.0, min(1.0, cos_val))
    return math.degrees(math.acos(cos_val))

def scout_branch(source_oid, target_oid, river_data, downstream_targets, max_depth=5):
    score = 0
    current = target_oid
    src_norm = river_data[source_oid]["norm"]
    src_fclass = river_data[source_oid]["fclass"]
    
    dev = calc_deviation(river_data[source_oid], river_data[target_oid])
    if dev > 30:
        score -= (dev * 2) 
    
    for _ in range(max_depth):
        if current not in river_data: break 
        
        curr_data = river_data[current]
        score += curr_data["len"]
        
        if src_norm and curr_data["norm"] == src_norm:
            score += 10000  
            
        if src_fclass and curr_data["fclass"] == src_fclass:
            score += 5000
            
        nxt_opts = [t for t in downstream_targets.get(current, []) if t in river_data]
        if not nxt_opts: 
            break 
            
        current = max(nxt_opts, key=lambda t: (
            src_norm != "" and river_data[t]["norm"] == src_norm,
            src_fclass != "" and river_data[t]["fclass"] == src_fclass,
            river_data[t]["len"]
        ))
            
    return score

def get_param(index, default_value):
    try:
        val = arcpy.GetParameterAsText(index)
        if not val: return default_value
        if isinstance(default_value, float): return float(val.replace(',', '.'))
        if isinstance(default_value, int): return int(float(val.replace(',', '.')))
        return val
    except:
        return default_value

# ==============================================================
# PARAMETRY A PROSTŘEDÍ
# ==============================================================
SRC_LINES       = arcpy.GetParameterAsText(0) or r"C:\Users\adams\Documents\BP_project_GIS\DATA IN\gis_osm_waterways_free_1.shp"
OUT_FOLDER      = arcpy.GetParameterAsText(1) or r"C:\Users\adams\Documents\BP_project_GIS\DATA"
TARGET_SCALE    = get_param(2, 100000.0)

MAIN_RIVER_NAME = arcpy.GetParameterAsText(4) or "Morava"
main_river_words = get_name_words(MAIN_RIVER_NAME.replace(',', ' '))

cutoff_raw      = arcpy.GetParameterAsText(6) or "umělý kanál, přivaděč, náhon, kanál"
CUTOFF_NAMES    = [w for n in cutoff_raw.split(',') for w in get_name_words(n)] if cutoff_raw else ["kanal", "privadec", "nahon"]
SUSPICIOUS_LIMIT = 50000.0 

TARGET_WKID     = 5514
NAME_FIELD      = "name"

arcpy.env.overwriteOutput = True
arcpy.env.parallelProcessingFactor = "100%"

if not os.path.exists(OUT_FOLDER):
    os.makedirs(OUT_FOLDER)

ts = datetime.now().strftime("%Y%m%d_%H%M%S")
safe_river = list(main_river_words)[0].capitalize() if main_river_words else "Reka"
gdb_path = os.path.join(OUT_FOLDER, f"Povodi_{safe_river}_{ts}.gdb")
arcpy.management.CreateFileGDB(OUT_FOLDER, f"Povodi_{safe_river}_{ts}.gdb")
arcpy.env.workspace = gdb_path

log(f"=== SPUŠTĚNO PRO MĚŘÍTKO 1:{int(TARGET_SCALE)} (Povodí: {MAIN_RIVER_NAME}) ===")

# =========================================================
# FÁZE A: PERFEKTNÍ SEGMENTACE (SOUTOKY)
# =========================================================
log("\n[FÁZE A] Projekce a planární rozbití (Feature To Line)...")
lines_proj = os.path.join(gdb_path, "L01_Projected")
arcpy.management.Project(SRC_LINES, lines_proj, arcpy.SpatialReference(TARGET_WKID))

lines_fc = os.path.join(gdb_path, "L02_Planarized")
arcpy.management.FeatureToLine(lines_proj, lines_fc, attributes="ATTRIBUTES")

# =========================================================
# FÁZE B: STRIKTNĚ ORIENTOVANÝ GRAF A "SCOUT"
# =========================================================
log(f"\n[FÁZE B] 1) Extrakce geometrie a budování jednosměrného grafu (Konec -> Začátek)...")
temp_end_pts = os.path.join(gdb_path, "temp_end_pts")
arcpy.management.FeatureVerticesToPoints(lines_fc, temp_end_pts, "END")

temp_start_pts = os.path.join(gdb_path, "temp_start_pts")
arcpy.management.FeatureVerticesToPoints(lines_fc, temp_start_pts, "START")

temp_near = os.path.join(gdb_path, "temp_near_directed")
arcpy.analysis.GenerateNearTable(temp_end_pts, temp_start_pts, temp_near, "0.1 Meters", closest="ALL")

endpt_to_seg = {pt: orig for pt, orig in arcpy.da.SearchCursor(temp_end_pts, ["OID@", "ORIG_FID"])}
startpt_to_seg = {pt: orig for pt, orig in arcpy.da.SearchCursor(temp_start_pts, ["OID@", "ORIG_FID"])}

river_data = {}
valid_class_fields = ['fclass', 'waterway', 'type', 'class']
class_field = next((f.name for f in arcpy.ListFields(lines_fc) if f.name.lower() in valid_class_fields), None)

search_fields = ["OID@", NAME_FIELD, "SHAPE@LENGTH", "SHAPE@"]
if class_field:
    search_fields.insert(2, class_field)

with arcpy.da.SearchCursor(lines_fc, search_fields) as cur:
    for row in cur:
        oid, nm = row[0], row[1]
        fcls = str(row[2]).lower() if class_field and row[2] else ""
        length = row[3] if class_field else row[2]
        geom = row[4] if class_field else row[3]
        first_pt = geom.firstPoint if geom else None
        last_pt = geom.lastPoint if geom else None
        
        river_data[oid] = {
            "name": nm if nm else "", "norm": normalize_name(nm), "words": get_name_words(nm), 
            "fclass": fcls, "len": length,
            "first": (first_pt.X, first_pt.Y) if first_pt else (0,0),
            "last": (last_pt.X, last_pt.Y) if last_pt else (0,0)
        }

downstream_targets = defaultdict(list)
upstream_sources = defaultdict(list)

with arcpy.da.SearchCursor(temp_near, ["IN_FID", "NEAR_FID"]) as cur:
    for pt_end, pt_start in cur:
        source_seg = endpt_to_seg.get(pt_end)
        target_seg = startpt_to_seg.get(pt_start)
        if source_seg and target_seg and source_seg != target_seg:
            if source_seg in river_data and target_seg in river_data:
                downstream_targets[source_seg].append(target_seg)
                upstream_sources[target_seg].append(source_seg)

primary_downstream = {}
primary_upstream = defaultdict(list)

log("[FÁZE B] 2) Výběr hlavních ramen (Scout: Jméno + Typ toku + Délka)...")
for source, targets in downstream_targets.items():
    valid_targets = [t for t in targets if t in river_data and river_data[t]["len"] > 5.0]
    if not valid_targets: valid_targets = [t for t in targets if t in river_data]
    if not valid_targets: continue
        
    if len(valid_targets) == 1:
        best_target = valid_targets[0]
    else:
        scores = {t: scout_branch(source, t, river_data, downstream_targets) for t in valid_targets}
        best_target = max(scores, key=scores.get)
    
    primary_downstream[source] = best_target
    primary_upstream[best_target].append(source)

log(f"[FÁZE B] 3) Izolace povodí '{MAIN_RIVER_NAME}' (Tvoje původní stabilní metoda)...")
main_segs = [oid for oid, d in river_data.items() if d["words"].intersection(main_river_words)]
if not main_segs: raise ValueError(f"Řeka {MAIN_RIVER_NAME} nenalezena v datech!")

main_segs_set = set(main_segs)
start_node = max(main_segs, key=lambda oid: river_data[oid]["len"])

curr = start_node
visited = set([curr])
while True:
    nxt = primary_downstream.get(curr)
    if not nxt or nxt not in main_segs_set or nxt in visited:
        break
    visited.add(nxt)
    curr = nxt
    
outlet_oid = curr
stack = [outlet_oid]
basin_nodes = set([outlet_oid])

while stack:
    curr = stack.pop()
    for c in upstream_sources[curr]:
        if c not in basin_nodes:
            basin_nodes.add(c)
            stack.append(c)

# =========================================================
# FÁZE C: TOPOLOGICKÉ ŘAZENÍ A AKUMULACE
# =========================================================
log("\n[FÁZE C] 1) Rychlé Kahnův řazení (O(1) Cycle Breaker)...")

for field in ["upstr_len", "render_val", "filter_val"]:
    arcpy.management.AddField(lines_fc, field, "DOUBLE")
for field in ["shreve_val", "strahler_val"]:
    arcpy.management.AddField(lines_fc, field, "LONG")

in_degree = {n: 0 for n in basin_nodes}
for n in basin_nodes:
    for c in primary_upstream[n]:
        if c in basin_nodes: in_degree[n] += 1

queue = deque([n for n in basin_nodes if in_degree[n] == 0])
topo_order = []
remaining = set(basin_nodes)

while remaining:
    if not queue:
        forced_node = next(iter(remaining))
        in_degree[forced_node] = 0
        queue.append(forced_node)
        
    curr = queue.popleft()
    if curr not in remaining: continue
        
    remaining.remove(curr)
    topo_order.append(curr)
        
    target = primary_downstream.get(curr)
    if target and target in remaining:
        in_degree[target] -= 1
        if in_degree[target] == 0:
            queue.append(target)

log("[FÁZE C] 2) Matematická akumulace v kostře (Domino efekt)...")
upstr_len, shreve_vals, strahler_vals = {}, {}, {}
oid_syms, oid_filters = {}, {}

# Agresivní prahy pro kartografický filtr
LEN_THRESHOLD = TARGET_SCALE * 0.08
SHR_THRESHOLD = max(1.0, TARGET_SCALE / 5000.0)

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

    l_val, sh_val, st_val = upstr_len[node], shreve_vals[node], strahler_vals[node]
    
    # NOVÝ VZOREC: sqrt(Shreve * Upstream Length)
    render_val = math.sqrt(max(1.0, sh_val) * max(1.0, l_val))
    
    # Filtrace: 1.0 (Ponechat) / 0.0 (Smazat)
    is_main = bool(river_data[node]["words"].intersection(main_river_words))
    f_val = 1.0 if (is_main or l_val >= LEN_THRESHOLD or sh_val >= SHR_THRESHOLD) else 0.0
    
    is_canal = any(c in river_data[node]["words"] for c in CUTOFF_NAMES)
    if is_canal and l_val > SUSPICIOUS_LIMIT:
        render_val, f_val = 0.0, 0.0

    oid_syms[node], oid_filters[node] = render_val, f_val

log("[FÁZE C] 3) Vizuální záchrana ostrovů (Nová multiplikační Ink Propagation)...")
ink_queue = deque()

for node in topo_order:
    alfa_target = primary_downstream.get(node)
    all_targets = [t for t in downstream_targets.get(node, []) if t in basin_nodes]
    
    for t in all_targets:
        if t != alfa_target:
            # Vedlejší rameno dostane 15 % hodnoty vizuální mohutnosti uzlu nad ním
            side_multiplier = 0.15 / max(1, len(all_targets) - 1)
            pass_sym = oid_syms.get(node, 0) * side_multiplier
            pass_flt = oid_filters.get(node, 0)
            ink_queue.append((t, pass_sym, pass_flt))

while ink_queue:
    curr, sym, flt = ink_queue.popleft()
    
    if sym > oid_syms.get(curr, -1):
        oid_syms[curr] = sym
        oid_filters[curr] = max(flt, oid_filters.get(curr, 0))
        
        for t in downstream_targets.get(curr, []):
            if t in basin_nodes:
                ink_queue.append((t, sym, flt))

log("  -> Zápis dat do tabulky atributů...")
with arcpy.da.UpdateCursor(lines_fc, ["OID@", "upstr_len", "shreve_val", "strahler_val", "render_val", "filter_val"]) as cur:
    for row in cur:
        oid = row[0]
        if oid in basin_nodes:
            row[1], row[2], row[3] = upstr_len.get(oid, 0), shreve_vals.get(oid, 0), strahler_vals.get(oid, 0)
            row[4], row[5] = oid_syms.get(oid, 0), oid_filters.get(oid, 0)
            cur.updateRow(row)
        else:
            cur.deleteRow()

# =========================================================
# FÁZE D: DYNAMICKÝ VÝBĚR (AGRESIVNÍ OŘEZ)
# =========================================================
log(f"\n[FÁZE D] Ořez sítě pro měřítko 1:{int(TARGET_SCALE)}...")
final_carto_fc = os.path.join(gdb_path, f"Povodi_{safe_river}_{int(TARGET_SCALE)}")

arcpy.analysis.Select(lines_fc, final_carto_fc, "filter_val >= 0.5")

try: 
    for tmp in [temp_end_pts, temp_start_pts, temp_near, lines_proj, lines_fc]:
        arcpy.management.Delete(tmp)
except: pass

# =========================================================
# FÁZE E: AUTOMATICKÁ KARTOGRAFIE 
# =========================================================
try:
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    m = aprx.activeMap
    
    if m:
        log("\n[FÁZE E] Vykreslování sítě (Pseudo-Proportional 0.1 - 7.0 pt)...")
        layer = m.addDataFromPath(final_carto_fc)
        sym = layer.symbology
        
        if hasattr(sym, 'updateRenderer'):
            sym.updateRenderer('GraduatedSymbolsRenderer')
            sym.renderer.classificationField = 'render_val'
            
            # API Trik: 100 tříd dělá vizuálně dokonalý Proportional efekt
            sym.renderer.breakCount = 100 
            
            sym.renderer.minimumSymbolSize = 0.1 
            sym.renderer.maximumSymbolSize = 7.0 
            
            color = {'RGB': [0, 92, 230, 100]}
            for brk in sym.renderer.classBreaks:
                try: brk.symbol.color = color
                except: pass
                
            layer.symbology = sym
            
        log(f"✅ Vrstva {layer.name} je připravena v mapě!")

except Exception as e:
    log(f"  [Info] Symbologii nastav ručně. (Běží mimo GUI: {e})")

log("\n" + "="*70)
log("ně d VÝPOČET KOMPLETNĚ DOKONČEN!")
log("="*70)