rm(list = ls(all=TRUE))

if (!require(tidyverse)) install.packages("tidyverse")
# if (!require(Seurat)) install.packages("Seurat")
if (!require(plotly)) install.packages("plotly")

library(tidyverse)
# library(Seurat)
library(plotly)

final_data <- tibble(
  X1 = numeric(),
  X2 = numeric(),
  epoch = numeric()
)

for ( i in list.files('../../../PoincareMaps/results/', pattern = 'seed0_epoch=[0-9]{2,4}.csv', full.name = TRUE) ) {
  # i = "../../../PoincareMaps/results//TandRepInfoShort_PM15sigma=1.00gamma=2.00minkowskipca=20_seed0_epoch=1000.csv"
  epoch <- gsub(pattern = "../../../PoincareMaps/results//TandRepInfoShort_PM15sigma=1.00gamma=2.00minkowskipca=20_seed0_epoch=", replacement = '', basename(i), fixed = TRUE) %>%
    gsub(pattern = 'TandRepInfoShort_PM15sigma=1.00gamma=2.00minkowskipca=20_seed0_epoch=', replacement = '') %>%
    gsub(pattern = '\\.csv', replacement = '') %>%
    as.integer()
  temp_data <- read_csv(i, col_names = FALSE, col_types = cols()) %>%
    mutate(
      epoch = epoch,
    )
  final_data <- bind_rows(final_data, temp_data)
}

final_data <- final_data %>% arrange(epoch)

glimpse(final_data)

to_plot = final_data[final_data$epoch %in% seq(10, 5000, 100),]

final_data$label = 1

fig <- plot_ly(
  to_plot,
  type = 'scatter',
  mode = 'markers',
  x = ~X1,
  y = ~X2,
  # color = ~label,
  frame = ~epoch,
  marker = list(size = 8),
  hoverinfo = 'text',
  text = ''
) %>%
  add_markers(
    data = to_plot %>%
      dplyr::select(epoch) %>%
      distinct() %>%
      mutate(x = 0, y = 0, label = 'center'),
    x = ~x,
    y = ~y,
    frame = ~epoch,
    marker = list(
      size = 10,
      line = list(
        color = 'black',
        width = 2
      ),
      symbol = 'x-thin'
    ),
    hoverinfo = 'skip',
    text = ~label,
    showlegend = FALSE
  ) %>%
  layout(
    title = 'Dynamics of training a Poincaré map',
    shapes = list(
      list(
        type = 'circle',
        xref = 'x',
        x0 = -1,
        x1 = 1,
        yref = 'y',
        y0 = -1,
        y1 = 1,
        fillcolor = NA,
        line = list(color = 'black'),
        opacity = 1
      )
    ),
    xaxis = list(
      range = c(-1, 1),
      showgrid = FALSE,
      showticklabels = FALSE,
      zeroline = FALSE,
      title = ''
    ),
    yaxis = list(
      range = c(-1, 1),
      showgrid = FALSE,
      showticklabels = FALSE,
      zeroline = FALSE,
      title = '',
      scaleanchor = 'x'
    )
  )

htmlwidgets::saveWidget(fig, 'poincare_map.html')
