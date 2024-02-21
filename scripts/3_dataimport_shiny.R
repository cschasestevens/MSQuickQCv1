#### Data Import ####

## format input data

names(list.data) <- list.qc.data[["File.ID"]]



## Assign df numbers and count sample numbers

list.qc.nos <- seq(nrow(list.qc.data))

list.qc.nos2 <- data.frame(Sample.No = t(as.data.frame(lapply(list.qc.nos,
                                                                 function(x) fun.count.samp(list.data[[x]]$Datasheet,
                                                                                            list.qc.data,
                                                                                            x)
                                                             )
                                                      )
                                        )
                          )

list.qc.data <- cbind(list.qc.data,
                      list.qc.nos2)





## Data imputation


### Replace NA and zero values with 1/10 of the lowest non-zero value for log transformation

list.data[["Imputed"]] <- lapply(list.qc.nos,
                                    function(x) 
                                      as.data.frame(
                                        lapply(
                                          list.data[[x]][["Datasheet"]],
                                          function(y) {

                                            ifelse(y == 0 |
                                                     is.na(y),
                                                   0.1*min(x[x > 0]),
                                                   y)
                                            }
                                          )
                                        )
                                 )

names(list.data[["Imputed"]]) <- list.qc.data[["File.ID"]]

















