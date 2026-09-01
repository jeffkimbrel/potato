# Simple plot of v2 graph

Simple plot of v2 graph

## Usage

``` r
plot_potato(potato_or_graph, interactive = TRUE, layout = "fr")
```

## Arguments

- potato_or_graph:

  Either a potato_v2 object, path to v2 JSON, or igraph from
  build_graph()

- interactive:

  Use interactive visNetwork plot (TRUE) or static ggraph (FALSE)

- layout:

  Layout algorithm for static plot (e.g., "fr", "kk", "circle", "tree",
  "grid"). Default: "fr"
