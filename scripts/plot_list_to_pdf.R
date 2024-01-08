# Function to save a list of ggplot files to PDF
plot_list_to_pdf <- function(plot_list, output_path, width = 21, height = 28, units = 'cm') {
  if (is.null(output_path)) {stop('Please provide an output path for the PDF file.')}

  # Convert sizes if needed
  if (units == 'cm') {
    width = width / 2.54
    height = height / 2.54
  }

  # Save plots to PDF
  pdf(file = output_path, width = width, height = height)
  if (class(plot_list)[1] == 'gg') {
    print(plot_list)
  } else {
    for (i in seq_along(plot_list)) {
      print(plot_list[[i]])
    }
  }
  dev.off()
}