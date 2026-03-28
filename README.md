# OSM River Network Generalizer (OSMgen)

An ArcGIS Pro geoprocessing tool developed as part of a Bachelor's thesis (2026). It takes raw, topologically inconsistent OpenStreetMap (OSM) waterway data and transforms it into a hydrologically correct, directed acyclic graph (DAG). The final output consists of generalized 2D river polygons (tapered based on stream order) and 1D centerlines, optimized for small-scale mapping (e.g., 1:100,000).

[🇨🇿 Přejít na českou verzi](#česká-verze)

## The Problem
Crowdsourced VGI data (like OSM) is messy. It contains missing topological nodes, artificial canals that break the main river trunk, unconnected river islands, and microscopic gaps at confluences. Standard geometric cutoffs and basic GIS tracing tools usually fail or destroy the actual network structure.

## Under the Hood (Key Algorithms)
This tool abandons simple geometric tracing in favor of graph theory and heuristics:
* **Tuple-Based Scout (Pathfinding):** At every bifurcation, the algorithm looks ahead up to 100 segments to evaluate the best path. It uses tuple scoring (exact name match > partial name match > waterway type > length) to prevent the main river from being diverted into a random artificial canal.
* **1.5m Topological Magnet:** Replaces exact spatial intersection. It bridges micro-gaps caused by careless manual digitizing or coordinate system projections.
* **Kahn’s Algorithm & Shreve Magnitude:** Topologically sorts the network from sources to the outlet to calculate the Shreve stream order.
* **Reverse Propagation (Backbone Protection):** Propagates the maximum Shreve magnitude from the river mouth back up the main trunk to the source. This prevents the headwaters of major rivers from being deleted during length-based truncation.
* **BFS Ink Propagation:** Uses a queue (FIFO) to push a percentage of "visual weight" into secondary channels, preventing river islands and anastomosing structures from being wiped out by filters.
* **River Tapering & Generalization:** Uses Douglas-Peucker and PAEK for cartographic smoothing. The polygon width is dynamically generated using $W = k \cdot \sqrt{Shreve}$ via a multi-threaded Pairwise Buffer.

## Usage
1. Download `OSMgen.atbx` and `OSMgen.py`.
2. In ArcGIS Pro, right-click **Toolboxes** -> **Add Toolbox** and select the `.atbx` file.
3. Open the tool, input your raw OSM lines, specify the target scale (e.g., 100000), target river name, and the output spatial reference. 

---

<a name="česká-verze"></a>
# OSM River Network Generalizer (OSMgen) 🇨🇿

Geoprocessingový nástroj pro ArcGIS Pro vyvinutý v rámci bakalářské práce (2026). Nástroj zpracovává surová a topologicky nekonzistentní data vodních toků z OpenStreetMap (OSM) a převádí je na orientovaný acyklický graf (DAG). Výstupem jsou generalizované 2D polygony říčních koryt (s proměnnou šířkou) a 1D osy optimalizované pro malá mapová měřítka (např. 1:100 000).

## Řešený problém
Crowdsourcovaná data (VGI) z OSM obsahují množství chyb: chybějící uzly na křižovatkách, umělé kanály přerušující hlavní tok, oddělená ramena ostrovů a mikroskopické mezery. Standardní plošné ořezy a jednoduché trasovací GIS nástroje na těchto datech většinou selhávají nebo ničí celistvost povodí.

## Jak to funguje pod kapotou
Skript opouští jednoduché geometrické trasování a spoléhá na teorii grafů a heuristiku:
* **Tuple-Based Scout (Pathfinding):** Na každé bifurkaci se algoritmus podívá až 100 segmentů dopředu. Trasování hodnotí přes n-tice (přesná shoda jména > částečná shoda > typ toku > délka), aby hlavní řeka neuhnula do bezejmenného náhonu.
* **Topologický magnet (1,5 m):** Namísto exaktního geometrického průniku skript spojuje uzly s mírnou tolerancí. Přemosťuje tak mikromezery vzniklé špatnou digitalizací nebo projekčními posuny.
* **Kahnův algoritmus a Shreveho řád:** Topologicky seřadí síť od pramenů k ústí pro matematickou akumulaci magnitudy.
* **Zpětná propagace (Ochrana páteře):** Skript pošle maximální Shreveho řád od ústí zpět proti proudu až k prameni hlavní řeky. To zabraňuje smazání pramenných úseků při plošném ořezu krátkých linií.
* **BFS Ink Propagation:** Pomocí fronty (FIFO) propouští část "vizuální váhy" do vedlejších ramen. Zabraňuje tak smazání říčních ostrovů a anastomózních struktur.
* **Generalizace a tvorba koryta:** Skript síť vyhladí (Douglas-Peucker + PAEK). Šířka finálního polygonu je počítána dynamicky vztahem $W = k \cdot \sqrt{Shreve}$ za pomoci rychlého nástroje Pairwise Buffer.

## Použití
1. Stáhněte si `OSMgen.atbx` a `OSMgen.py`.
2. V ArcGIS Pro klikněte pravým tlačítkem na **Toolboxes** -> **Add Toolbox** a připojte stažený `.atbx` soubor.
3. Spusťte nástroj, vložte surové OSM linie, zadejte cílové měřítko (např. 100000), název zájmové řeky a cílový souřadnicový systém.

---
**Autor:** Adam Sojka 
**Licence:** Nástroj je k dispozici open-source pro analytické a akademické účely.
