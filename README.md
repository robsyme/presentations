# presentations

Slide decks and small interactive tools, published at
<https://robsyme.github.io/presentations/>.

Each one is a directory of plain static files, with an `index.html` and its own
assets. They do not share code, so one can be changed or deleted without touching the
others. A deck holds the slide markup and the small runtime that renders it, and
nothing in it is built or compiled. A tool may be built from source kept elsewhere;
its directory holds the built page.

## Presentations

| Date | Deck |
|------|------|
| September 2026 | [Nextflow: the ideas before the syntax](2026-09-nextflow-ideas-before-syntax/) |

## Tools

| Date | Tool |
|------|------|
| November 2026 | [Barcode Contact Map](2026-11-barcode-contact-map/): viewer for linked-read barcode co-occurrence; [guide](2026-11-barcode-contact-map/guide.html) and `preprocess.py` alongside |
| October 2026 | [Québec knowledge test prep](2026-10-quebec-knowledge-test/): interactive study guide |

## Running one locally

The decks use `fetch()` to pull in their components, so they need to be served over
HTTP rather than opened from the filesystem:

```
python3 -m http.server 8000
```

Then open <http://localhost:8000/>.

## Presenting

Right arrow, down arrow, space or PageDown advance a slide; left, up or PageUp go
back. `Home` and `End` jump to the first and last slide, `R` restarts, and the number
keys jump to slides 1–9 (`0` for 10). Slides carrying a Rive animation advance one
step per click.

Speaker notes live in the `data-speaker-notes` attribute of each `<section>` in the
deck's `index.html`. Adding `#7` to a deck's URL opens it at slide 7, and the
browser's Print → Save as PDF gives one slide per page.

Every runtime dependency is vendored under each deck's `vendor/` directory, so a deck
will present on a laptop with no network — you still need the local HTTP server above.

## Adding a deck

Copy the exported deck directory in under a `YYYY-MM-slug` name, rename its
`*.dc.html` to `index.html`, and add a row to the Presentations table above and an
entry under Presentations in the root `index.html`. Worth doing at the same time, since each export arrives fatter than it
needs to be:

- Drop fonts the deck does not use, and their `@font-face` rules.
- Drop background images and other assets nothing references.
- Vendor the CDN scripts (see `2026-09-nextflow-ideas-before-syntax/vendor/` and the
  `window.__resources` block at the top of its `index.html`) if the deck needs to
  survive a bad conference network.
- Delete the design system's `_ds_bundle.js` and the `<script>` that loads it. It
  carries an older copy of `deck-stage.js` that wins the race to define the element,
  costing you the newer build's print-to-PDF and `#<n>` slide links.
- Turn off the editor's thumbnail rail. The deck in this repo does it with a short
  script that sets `no-rail` on `<deck-stage>` once the runtime has mounted it.

## Adding a tool

Put the built page in a `YYYY-MM-slug` directory as `index.html`, with any files it
links to (a guide, a downloadable script) beside it. Add a row to the Tools table above
and an entry under Tools in the root `index.html`. A tool should not load anything from
the network, so it keeps working offline like the decks.

`.nojekyll` at the repo root stops GitHub Pages from running Jekyll, which would
otherwise refuse to serve any directory whose name begins with an underscore.

## Fonts

The decks self-host [Inter](https://github.com/rsms/inter) and
[JetBrains Mono](https://github.com/JetBrains/JetBrainsMono), both under the SIL Open
Font License. Licensed fonts are not redistributed here; where a deck's design system
names one, the text falls through to the system stack.
