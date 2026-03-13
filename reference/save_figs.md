# Save Figures

Save Figures

## Usage

``` r
save_figs(
  name,
  fig,
  width = 6,
  height = 6,
  plot_dir = file.path(here::here(), "analysis/plots"),
  pdf_plot = TRUE,
  font_family = "Helvetica",
  res = 300,
  ...
)
```

## Arguments

- name:

  Name of figure

- fig:

  ggplot or similar figure object

- width:

  Width of plot in inches. Default = 6

- height:

  Height of plot in inches. Default = 6

- plot_dir:

  Plotting directory. Defaults to "analysis/plots"

- pdf_plot:

  Logical for plotting pdf too. Default = TRUE

- font_family:

  If specified, sets all font family. Default = NULL

- res:

  Image resolution in dpi. Default = 300

- ...:

  Other parameters to pass to ragg::agg_png
