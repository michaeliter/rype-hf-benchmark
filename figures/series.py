"""Single registry of which approaches the figures show, and how.

To add a tool: get its rows into data/<mod>/ (eval/eval_op.py), give it an
entry in COLOR/LABEL/MARKER below, and add its key to EXTERNAL. Every figure
picks it up from here; plot order follows EXTERNAL.
"""

RYPE = ["rype", "rype-mm2"]
EXTERNAL = ["cleanifier", "cleanifier-filter", "movi2-alone", "mm2-mm2-movi"]

# Paul Tol palette; cleanifier-filter is Tol purple, chosen to pass CVD
# separation against its neighbors.
COLOR = {"rype": "#0077BB", "rype-mm2": "#009988", "mm2-mm2-movi": "#CC3311",
         "movi2-alone": "#EE7733", "cleanifier": "#BBBBBB",
         "cleanifier-filter": "#AA3377"}
LABEL = {"rype": "RYpe", "rype-mm2": "RYpe + mm2", "mm2-mm2-movi": "mm2→mm2→Movi (SOTA)",
         "movi2-alone": "Movi 2 alone (no mm2)",
         "cleanifier": "cleanifier (exact)", "cleanifier-filter": "cleanifier (probabilistic)"}
MARKER = {"rype": "o", "rype-mm2": "s", "mm2-mm2-movi": "^", "movi2-alone": "D",
          "cleanifier": "P", "cleanifier-filter": "X"}
# Storage components (storage.tsv) for each window-independent tool's index.
INDEX_COMPONENTS = {
    "cleanifier": [("human-all97 exact (cleanifier)", "cleanifier_index")],
    "cleanifier-filter": [("human-all97 filter (cleanifier)", "cleanifier_filter_index")],
    "movi2-alone": [("hprc-yr1 (Movi 2)", "movi2_index")],
}

APPROACHES = RYPE + EXTERNAL
WINDOW_INDEPENDENT = set(EXTERNAL)
