"""Extract Parc <-> Bois-de-Boulogne train times from exo's train GTFS feed.

Usage: python3 build_timetable.py <unzipped-gtfs-dir> > timetable.json

The page does not read timetable.json. After a schedule change, copy one weekday
service's outbound/inbound times into the TT constant in index.html as
[departs, arrives, train] triples (minutes after midnight).

Feed: https://exo.quebec/uploads/ressources-telechargeables/a-propos/donnees-ouvertes/trains/google_transit.zip
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ORIGIN = "Gare Parc"
DEST = "Gare Bois-de-Boulogne"


def rows(gtfs, name):
    with open(gtfs / name, newline="", encoding="utf-8-sig") as f:
        yield from csv.DictReader(f)


def hhmm_to_min(t):
    h, m, s = (int(x) for x in t.strip().split(":"))
    return h * 60 + m + round(s / 60)


def main(gtfs):
    stops = {r["stop_id"]: r["stop_name"].strip() for r in rows(gtfs, "stops.txt")}
    weekday_services = {
        r["service_id"]: {"start": r["start_date"], "end": r["end_date"]}
        for r in rows(gtfs, "calendar.txt")
        if all(r[d] == "1" for d in ("monday", "tuesday", "wednesday", "thursday", "friday"))
    }
    trips = {
        r["trip_id"]: r
        for r in rows(gtfs, "trips.txt")
        if r["service_id"] in weekday_services
    }

    times = defaultdict(dict)  # trip_id -> stop_name -> (arr, dep, pickup_ok, dropoff_ok)
    for r in rows(gtfs, "stop_times.txt"):
        if r["trip_id"] not in trips:
            continue
        name = stops[r["stop_id"]]
        if name not in (ORIGIN, DEST):
            continue
        times[r["trip_id"]][name] = {
            "arr": hhmm_to_min(r["arrival_time"]),
            "dep": hhmm_to_min(r["departure_time"]),
            "pickup": r["pickup_type"] != "1",
            "dropoff": r["drop_off_type"] != "1",
        }

    out = {
        sid: {"start": v["start"], "end": v["end"], "outbound": [], "inbound": []}
        for sid, v in weekday_services.items()
    }
    for tid, st in times.items():
        if ORIGIN not in st or DEST not in st:
            continue
        trip = trips[tid]
        o, d = st[ORIGIN], st[DEST]
        entry = {"train": trip["trip_short_name"], "headsign": trip["trip_headsign"]}
        if o["dep"] < d["arr"] and o["pickup"] and d["dropoff"]:
            out[trip["service_id"]]["outbound"].append({**entry, "dep": o["dep"], "arr": d["arr"]})
        elif d["dep"] < o["arr"] and d["pickup"] and o["dropoff"]:
            out[trip["service_id"]]["inbound"].append({**entry, "dep": d["dep"], "arr": o["arr"]})

    for svc in out.values():
        svc["outbound"].sort(key=lambda x: x["dep"])
        svc["inbound"].sort(key=lambda x: x["dep"])
    json.dump(out, sys.stdout, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
