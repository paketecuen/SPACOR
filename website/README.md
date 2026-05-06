# SPACOR Website

Static public website for `spacor.ual.es`.

The site is intentionally static: it does not need a database, PHP framework, or
private datasets.  Deploy the contents of this folder to the web root served by
the server, currently expected to be:

```text
/data/spacor.ual.es/web/public/
```

The page structure is:

- `index.html`: public scientific landing page.
- `extended.html`: online supplementary material and extended figures.
- `assets/css/site.css`: site styles.
- `assets/figures/`: selected public figures generated from synthetic campaigns.

Generated or private data should not be copied into this folder.
