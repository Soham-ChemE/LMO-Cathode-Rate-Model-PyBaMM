# QA for the live page (2026-11-04)

Served from docs/ with python -m http.server, loaded in headless Chromium (software WebGL, so render timing is not representative).

## Lighthouse (desktop preset)

Performance 98, Accessibility 98, Best Practices 100, SEO 100

| Metric | Value |
|---|---|
| First Contentful Paint | 0.9 s |
| Largest Contentful Paint | 0.9 s |
| Speed Index | 1.1 s |
| Total Blocking Time | 0 ms |
| Cumulative Layout Shift | 0 |
| Time to Interactive | 0.9 s |
| Avoids enormous network payloads | Total size was 980 KiB |

Report: lighthouse.report.html, lighthouse.report.json.

## Headless-browser check

| View | Load (ms) | Transfer (bytes) | Images loaded | Canvases | Failed requests | Console errors |
|---|---|---|---|---|---|---|
| desktop | 1335 | 2,873,762 | 4/4 | heroCanvas 1440x900 | 0 | 0 |
| mobile | 1252 | 2,873,762 | 4/4 | heroCanvas 375x812 | 0 | 0 |

Screenshots: 2026-11-04-desktop-hero.png, 2026-11-04-desktop-dark-ribbons.png, 2026-11-04-desktop-end.png, 2026-11-04-desktop-full.png, 2026-11-04-mobile-375-hero.png, 2026-11-04-mobile-375-full.png.
