# Function to rotate a data.frame and maintain names
rotate_df <- function(data){
  names_row <- rownames(data)
  names_col <- colnames(data)
  
  data_rotated <- data.frame(t(data))
  rownames(data_rotated) <- names_col
  colnames(data_rotated) <- names_row
  
  data_rotated
}