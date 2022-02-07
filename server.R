#NEW VERSION
server <- function(input, output) {

  #options(shiny.maxRequestSize=100*1024^2)
  # this prints a dot to the cosole every 10 seconds so that the screen doesn't time out
  autoInvalidate <- reactiveTimer(10000)
  observe({
    autoInvalidate()
    cat(".")
  })

  # tell the shinyhelper package what the file name of the help file is
  observe_helpers(help_dir = "HelpFile")

  observe({

    file1 = input$inFile
    if (is.null(file1)) {
      return(NULL)
    }

    data1 = read.csv(file1$datapath)
    #if(input$format == ".csv") data1 = read.csv(file1$datapath)
    #else data1 = read_excel(file1$datapath)

    if(length(unique(data1[,1]))!=length(data1[,1])) {

      n = as.character(data1[,1])
      dups = which(duplicated(n))
      output$warning = renderText("<b>These sample names are duplicated in your dataset:")
      output$warning2 = renderText(c("<B>", n[dups]))
      output$warning3 = renderText("<b>Please ensure that all sample names are unique. Refresh the program, and re-upload the file.")
    }

    else if ((nrow(data1) %% 2) != 0) {

      output$warning4 = renderText("<b>There is an odd number of samples in this dataset. If you have a set of replicate pairs, there should be an even number. Please correct this, refresh the program, and re-upload your file.")

    }

    else{
      # removes error warnings if there are any from a previous upload
      output$warning = renderText("")
      output$warning2 = renderText("")
      output$warning3 = renderText("")
      output$warning4 = renderText("")

      row.names(data1) <- data1[[1]]

      data1[,1] <- NULL

      #data1[,] <- lapply(data1[,], gsub, pattern = "$", replacement = "", fixed = TRUE)
      #data1[,] <- lapply(data1[,], gsub, pattern = "?", replacement = "1", fixed = TRUE)
      data1[,] <- sapply(data1[,], as.numeric)

      #data1[data1 > 1] = 1

      #round(data1,0)
      #observeEvent(input$preview, {output$table <- renderTable(rownames=T,{data1},digits = 0)

      output$downloadData <- downloadHandler(

        filename = function (){paste('ConsolidatedBinaryMatrix', '.csv', sep = '.')},
        content = function (file){write.csv(data1, file)}

      )

      #})

      answer = which(data1 != 0 & data1 != 1 & data1 != "?", arr.ind = T)


      observeEvent(input$check,{output$text1 <- renderTable(if(length(answer) > 0){answer} else {"None found"}, caption = "Unwanted values at:", caption.placement = getOption("xtable.caption.placement", "top"))})


      observeEvent(input$act, {

        # if(input$choice == 2)
        # {

        odd = data1[seq(1, nrow(data1), by = 2),]
        even = data1[seq(2, nrow(data1), by = 2),]
        new = odd+even

        new[,] <- lapply(new[,], gsub, pattern = "1", replacement = "?", fixed = TRUE)
        new[,] <- lapply(new[,], gsub, pattern = "2", replacement = "1", fixed = TRUE)
        nams = row.names(data1)
        samplenames = vector()
        for(i in 1:(length(nams)/2))
        { #divide by two, because the new matrix is half the size (due to the rep pairs being combined)
          samplenames[i] = paste(nams[i*2-1], nams[i*2], sep = '+')
        }

        row.names(new) = samplenames
        colnames(new) = colnames(data1)

        output$consod_done = renderText("<b>COMPLETE. READY FOR DOWNLOAD.")
        #output$table <- renderTable(rownames = T,{new})

        mismatch_err = matrix(nrow=nrow(new), ncol = 1)
        jacc_err = matrix(nrow=nrow(new), ncol = 1)


        for(i in 1:nrow(new)) {
          # find the number of 1s, Os and question marks
          ones = length(which(new[i,] == 1))
          zeroes = length(which(new[i,] == 0))
          questions = length(which(new[i,] == "?"))
          sum_bands = ones + questions

          mismatch_err[i,] = (questions/(questions + ones + zeroes))
          jacc_err[i,] = (questions/(questions + ones))

        }

        error_table = data.frame("Errors" = matrix(ncol = 2, nrow = 4))
        error_table[1,1] = "Average Euclidean Error:"
        error_table[1,2] = round(mean(mismatch_err[,1]),4)
        error_table[2,1] = "Euclidean error St. dev:"
        error_table[2,2] = round(sd(mismatch_err[,1]),4)
        error_table[3,1] = "Average Jaccard:"
        error_table[3,2] = round(mean(jacc_err[,1]),4)
        error_table[4,1] = "Jaccard error St.dev:"
        error_table[4,2] = round(sd(jacc_err[,1]),4)
        
        colnames(error_table) = c("Metric", "Value")

        observeEvent(input$repro,{

          output$text2 <- renderTable(error_table)

          #output$text2 <- renderTable(c("Average Euclidean Error:", round(mean(mismatch_err[,1]),4), "Average Jaccard:", round(mean(jacc_err[,1]),4), "Euclidean error St. dev:", round(sd(mismatch_err[,1]),4), "Jaccard error St.dev:", round(sd(jacc_err[,1]),4)), caption = "Error Rates:", caption.placement = getOption("xtable.caption.placement", "top"))

          output$download_errors <- downloadHandler(

            filename = function (){paste('Error_rates', 'csv', sep = '.')},
            content = function (file){write.csv(error_table, file, row.names = F)}
          )

        })

        observeEvent(input$summary, {


          nr_peaks = matrix(nrow = nrow(new), ncol = 1)


          for(i in 1:nrow(new)) {
            total = 0
            for(j in 1:ncol(new)) {if(new[i,j] == 1) total = total + 1}

            nr_peaks[i,] = total

          }

          summary_table = data.frame("Summary" = matrix(ncol = 2, nrow = 5))
          summary_table[1,1] = "Average no. peaks "
          summary_table[1,2] = round(mean(nr_peaks),4)
          summary_table[2,1] = "Standard deviation"
          summary_table[2,2] = round(sd(nr_peaks),4)
          summary_table[3,1] = "Max. no. peaks"
          summary_table[3,2] = max(nr_peaks)
          summary_table[4,1] = "Min. no. peaks"
          summary_table[4,2] = min(nr_peaks)
          summary_table[5,1] = "No. loci "
          summary_table[5,2] = ncol(new)
          
          colnames(summary_table) = c("Metric", "Value")

          output$text3 = renderTable(summary_table)
          #output$text3 = renderTable(c("Average no. peaks = ", round(mean(nr_peaks),4), "sd = ", round(sd(nr_peaks),4), "Max. no. peaks = ", max(nr_peaks), "Min. no. peaks = ", min(nr_peaks), "No. loci = ", ncol(new)))

          output$download_summary <- downloadHandler(

            filename = function (){paste('Error_rates', 'csv', sep = '.')},
            content = function (file){write.csv(summary_table, file, row.names = F)}
          )

        })

        observeEvent(input$remove, {

          percentage = input$jacc_remove
          new$jacc_error = jacc_err

          removed = subset(new, (new[,ncol(new)] > percentage))
          new = subset(new, (new[,ncol(new)] <= percentage))

          mismatch_err = matrix(data = NA, nrow=nrow(new), ncol = 1)
          jacc_err = matrix(data = NA, nrow=nrow(new), ncol = 1)

          for(i in 1:nrow(new)) {
            # find the number of 1s, Os and question marks
            ones = length(which(new[i,c(1:ncol(new))] == 1))
            zeroes = length(which(new[i,c(1:ncol(new))] == 0))
            questions = length(which(new[i,c(1:ncol(new))] == "?"))
            sum_bands = ones + questions
            mismatch_err[i,] = (questions/(questions + ones + zeroes))
            jacc_err[i,] = (questions/(questions + ones))

          }

          observeEvent(input$repro,{output$text2 <- renderTable(c("Average Euclidean Error:", round(mean(mismatch_err[,1]),4), "Average Jaccard:", round(mean(jacc_err[,1]),4), "Euclidean error St. dev:", round(sd(mismatch_err[,1]),4), "Jaccard error St.dev:", round(sd(jacc_err[,1]),4), "% of removed samples", round((nrow(removed)/(nrow(new) + nrow(removed))*100),4)), caption = "Error Rates:", caption.placement = getOption("xtable.caption.placement", "top"))})

          output$downloadData <- downloadHandler(

            filename = function (){paste('ConsolidatedBinaryMatrix', 'csv', sep = '.')},
            content = function (file){write.csv(new, file)}
          )
        })


        output$downloadData <- downloadHandler(

          filename = function (){paste('ConsolidatedBinaryMatrix', 'csv', sep = '.')},
          content = function (file){write.csv(new, file)}
        )


        #"draw tree" button
        # observeEvent(input$drawTree, {
        #
        #   new_mat = new
        #   # make the names shorter (here, only 10 characters long)
        #   row.names(new_mat) = row.names(new)
        #   newnames =substring(row.names(new_mat), 0, 50)
        #   row.names(new_mat) = newnames
        #   rownames(new_mat) = make.names(newnames, unique = T)
        #
        #   new_mat[new_mat=="?"] <- NA
        #   new_mat = as.data.frame(new_mat)
        #
        #   result <- pvclust(t(new_mat), method.dist = "binary", method.hclust="average", nboot=input$boot) # 'average' method is the UPGMA
        #
        #   output$upgmaTree = {renderPlot(plot(result, cex = 0.55, print.num = F, print.pv = T, cex.pv = 0.55))}
        #
        #
        #   output$downloadTree <- downloadHandler(
        #
        #     filename = function (){paste("UPGMAtree", "svg", sep = '.')},
        #
        #     content = function (file){svg(file)
        #
        #       plot(result, cex = 0.55, print.num = F, print.pv = T, cex.pv = 0.55)
        #       dev.off()
        #
        #     }
        #   )
        #
        # })



        # }  # end of if statement regarding data type


        # else {
        #   consod = data.frame()
        #
        #   for(i in 1:ncol(data1)){
        #
        #     avg = sum(data1[,i]/nrow(data1))
        #     if (avg > 0.5) consod[1,i] = paste("1")
        #     if (avg == 0.5) consod[1,i] = paste("?")
        #     if (avg < 0.5) consod[1,i] = paste("0")
        #   }
        #
        #   colnames(consod) = colnames(data1)
        #   row.names(consod) = "Combined replicates"
        #
        #   output$table <- renderTable(rownames = T,{consod})
        #
        #   # find the number of 1s, Os and question marks
        #   ones = length(which(consod == 1))
        #   zeroes = length(which(consod == 0))
        #   questions = length(which(consod == "?"))
        #   sum_bands = ones + questions
        #
        #   observeEvent(input$repro,{output$text2 <- renderTable(c("1s (%):", ones/sum_bands*100 , "? (%):", questions/sum_bands*100), caption = "Reproducibility:", caption.placement = getOption("xtable.caption.placement", "top"))})
        #
        #   output$downloadData <- downloadHandler(
        #
        #     filename = function (){paste('name', '.csv', sep = '.')},
        #     content = function (file){write.csv(consod, file)}
        #
        #   )
        #
        #   observeEvent(input$drawTree, {
        #     output$warning5 = renderText("<b> You can only create a clustering tree for data consisting of replicate pairs.")
        #
        #   })
        #
        # } # end of else statement

      })
    }
  })





  # for the separate uploading of data for an MDS plot:

  observe({

    file_mds = input$inFile_mds
    if (is.null(file_mds)) {
      return(NULL)
    }

    data_mds = read.csv(file_mds$datapath)

    if(length(unique(data_mds[,1]))!=length(data_mds[,1])) {

      n2 = as.character(data_mds[,1])
      dups2 = which(duplicated(n2))
      output$warn1 = renderText("<b>These sample names are duplicated in your dataset:")
      output$warn2 = renderText(c("<B>", n2[dups2]))
      output$warn3 = renderText("<b>Please ensure that all sample names are unique. Close/refresh the program, and re-upload the file.")
    }

    else{
      # clear warning messages if there were any
      output$warn1 = renderText("")
      output$warn2 = renderText("")
      output$warn3 = renderText("")


      row.names(data_mds) <- data_mds[[1]] # make the sample names rownames,
      data_mds[,1] <- NULL # and then remove the sample name column

      # make the names shorter (here, only 10 characters long)
      newnames =substring(row.names(data_mds), 0, 50)
      row.names(data_mds) = newnames
      data_mds[data_mds =="?"] <- NA


      data_mds = as.data.frame(data_mds)

      fac = factor(data_mds[[1]], levels = unique(data_mds[[1]])) # create a factor variable for groups

      #data_mds[[1]] = NULL # remove the group column from the dataframe

      colour_choices = c("dodgerblue", "black", "red", "green3", "orange", "darkblue", "gold2", "darkgreen", "darkred", "grey", "darkgrey", "magenta", "darkorchid", "purple", "brown", "coral3", "turquoise", "deeppink", "lawngreen", "darkred", "deepskyblue", "tomato", "yellow", "yellowgreen", "royalblue", "olivedrab", "midnightblue", "indianred1", "darkturquoise")
      #colour_choices = c(colors())
      shape_choices = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25")

      fact_table = data.frame(Groups = levels(fac), Colour_by_Group = colour_choices[1:length(levels(fac))], Point_Shape="16", stringsAsFactors=F)
      #Colour_by_Group = colors()[1:length(levels(fac))]
      #colour_choices = c(colors())


      output$grps = renderRHandsontable(
        rhandsontable(fact_table, height = "300")
        %>% hot_col(col = "Colour_by_Group", type = "dropdown", source = colour_choices)
        %>% hot_col(col = "Point_Shape", type = "dropdown", source = shape_choices)
      )
      
      if(input$display_labs == TRUE) mds_labels = row.names(data_mds)
      else mds_labels = ""
      
      pt_size = input$cexSize

      observeEvent(input$plotMDS, {

        dist_meth = input$dist_method

        d = dist((data_mds[,2:ncol(data_mds)]), method = dist_meth, diag = TRUE, upper = T)
        d = as.data.frame(as.matrix(d))
        d2 = as.dist(d)
        d2 = d2 + 0.01 # adding 0.01 here to cover for cases where there are identical sequences, leading to zero distances. Zero distances give the error "Warning: Error in isoMDS: zero or negative distance between objects x and y"

        k_val = input$k_value
        
        isoplot = isoMDS(d2, k = k_val)
        
        #isoplot_2 = data_mds[,2:ncol(data_mds)] %>% dist(method = dist_meth) %>% isoMDS(k = k_val) %>% .$points %>% tibble::as_tibble()
        #colnames(isoplot_2) = c("Dimension 1", "Dimension 2")
        
        # req(input$grps)
        colour_update = hot_to_r(input$grps)
        fact_table[,2] = colour_update[,2]
        fact_table[,3] = colour_update[,3]
        
        isoplot_df = as.data.frame(tibble::as_tibble( isoplot$points ))
        isoplot_df$colours = c(colour_update[[2]])[fac]
        isoplot_df$groups = fac
        
        isoplot_df$colours = factor(isoplot_df$colours, levels = unique(isoplot_df$colours))
        
        if (k_val == 2) colnames(isoplot_df) = c("Dimension 1", "Dimension 2", "colours", "groups")
        else colnames(isoplot_df) = c("Dimension 1", "Dimension 2", "Dimension 3", "colours", "groups")
        
        if(input$plot_components == "1 and 2"){
          x_comp = "Dimension 1"
          y_comp = "Dimension 2"
        }
        
        else if(input$plot_components == "1 and 3"){
          x_comp = "Dimension 1"
          y_comp = "Dimension 3"
        }
        
        else if(input$plot_components == "2 and 3"){
          x_comp = "Dimension 2"
          y_comp = "Dimension 3"
        }
       
        
        #colour_update contains the dataframe with 2 columns: group and colour

        output$mdsPlot = renderPlot({
          #par(mar = c(6.1, 8.1, 5.1, 7.1))
          
         pp = ggpubr::ggscatter(isoplot_df, 
                                x = x_comp, 
                                y = y_comp, 
                                label = mds_labels,
                                repel = input$repel_labs,
                                color = "groups",
                                palette = levels(isoplot_df$colours),
                                shape = c(as.numeric(colour_update[[3]]))[fac],
                                ellipse = input$display_ellipses,
                                ellipse.type = input$ellipse_type,
                                size = pt_size,
                                ellipse.alpha = input$ellipse_alpha,
                                star.plot = input$star_plot
                                )
         
         plot(pp)
         #eqscplot(isoplot$points, xlab = "Dimension 1", ylab = "Dimension 2", col = c(colour_update[[2]])[fac], pch = c(as.numeric(colour_update[[3]]))[fac], cex = pt_size)
          
        })
        
        output$download_dist_matrix = downloadHandler(
          filename = function (){paste('distance_matrix', 'csv', sep = '.')},
          content = function (file){write.csv(d, file)}
        )

        output$downloadMDS <- downloadHandler(

          filename = function (){paste("nMDS_Plot", "svg", sep = '.')},

          content = function (file){svg(file)
            pp = ggpubr::ggscatter(isoplot_df, 
                                   x = x_comp, 
                                   y = y_comp, 
                                   label = mds_labels,
                                   repel = input$repel_labs,
                                   color = "groups",
                                   palette = levels(isoplot_df$colours),
                                   shape = c(as.numeric(colour_update[[3]]))[fac],
                                   ellipse = input$display_ellipses,
                                   ellipse.type = input$ellipse_type,
                                   size = pt_size,
                                   ellipse.alpha = input$ellipse_alpha,
                                   star.plot = input$star_plot
            )
            plot(pp)
            dev.off()
          }
        )

        # this is within the plotMDS section, because it needs access to the isoplot object
        observeEvent(input$shep, {

          shep = Shepard(d2,isoplot$points, p=2)
          summ = summary(lm(shep$y~shep$x)) # R squared value
          r_sq = round(summ$r.squared, digits = 2)

          #grad = round(coef(lm(shep$y~shep$x))[2], digits = 2)
          #inrcpt = round(coef(lm(shep$y~shep$x))[1], digits = 2)

          output$shepPlot = renderPlot({
            plot(shep, pch=16, xlab = "Original data distances", ylab = "Ordination distances", main = c("R-squared = ", r_sq))
            abline(lm(shep$y~shep$x), col = "blue", lty = 2,lwd = 2)
          })

          output$download_shep <- downloadHandler(

            filename = function (){paste("Shepard_Plot", "svg", sep = '.')},

            content = function (file){svg(file)
              plot(shep, pch=16, xlab = "Original data distances", ylab = "Ordination distances")
              abline(lm(shep$y~shep$x), col = "blue", lty = 2,lwd = 2)
              dev.off()
            }
          )

        })

        observeEvent(input$scree, {

          scree.plot = function(d, k) {
            stresses=isoMDS(d, k=k)$stress
            for(i in rev(seq(k-1)))
              stresses=append(stresses,isoMDS(d, k=i)$stress)
            plot(seq(k),rev(stresses), type="b", xaxp=c(1,k, k-1), ylab="Stress (%)", xlab="Number of dimensions")

          }

          output$screePlot = renderPlot({
            scree.plot(d2, k=4)
            abline(b=0,h=15, col = "red", lty = 2)
          })

          output$download_scree <- downloadHandler(

            filename = function (){paste("Scree_Plot", "svg", sep = '.')},

            content = function (file){svg(file)
              scree.plot(d2, k=6)
              dev.off()
            }
          )

        })



      })

      observeEvent(input$filter_samples, {

        peak_thresh_remove = input$peak_thresh

        peak_record = matrix(data = NA, nrow(data_mds), ncol = 1)

        for(i in 1:nrow(data_mds)){
          count = length(which(data_mds[i,2:ncol(data_mds)]==1))
          peak_record[i,]=count
        }

        if(peak_thresh_remove > max(peak_record))
          output$msg = renderText(c("<b>Do not remove samples with peaks less than ", max(peak_record)+1, ", as that the maximum peak number detected in your data set is ", max(peak_record), "."))

        else {

          data_mds_tally = cbind(data_mds, peak_record)
          data_mds_keep = subset(data_mds_tally, peak_record >= peak_thresh_remove)
          data_mds_removed = subset(data_mds_tally, peak_record < peak_thresh_remove)
          data_mds_keep$peak_record = NULL
          data_mds_removed$peak_record = NULL
          names(data_mds_keep)[1]="Group"

          data_mds_keep[is.na(data_mds_keep)]<-"?"

          names(data_mds_removed)[1]="Group"

          output$download_filtered = downloadHandler(
            filename = function (){paste('subsetted_data', 'csv', sep = '.')},
            content = function (file){write.csv(data_mds_keep, file)}
          )

          output$download_removed = downloadHandler(
            filename = function (){paste('removed_data', 'csv', sep = '.')},
            content = function (file){write.csv(data_mds_removed, file)}
          )

          output$msg = renderText(c("<b>", nrow(data_mds_removed), "removed out of the ", nrow(data_mds), " samples provided."))

        }
      })


    } # end of else statement regarding duplicate row names


  })


  observe({

    file_upgma = input$inFile_upgma
    if (is.null(file_upgma)) {
      return(NULL)
    }


    observeEvent(input$drawTree, {

      data_upgma1 = read.csv(file_upgma$datapath)
      data_upgma2 = data_upgma1[,-1]
      row.names(data_upgma2) = data_upgma1[,1]
      new_names_upgma = substring(row.names(data_upgma2),0,50)
      row.names(data_upgma2) = new_names_upgma
      data_upgma2[data_upgma2=="?"] = NA
      data_upgma2 = as.data.frame(data_upgma2)
      result <- pvclust(t(data_upgma2), method.dist = "binary", method.hclust="average", nboot=input$boot) # 'average' method is the UPGMA

      output$upgmaTree = {renderPlot(plot(result, cex = 0.55, print.num = F, print.pv = T, cex.pv = 0.55))}

      output$downloadTree <- downloadHandler(
        filename = function (){paste("UPGMAtree", "svg", sep = '.')},
        content = function (file){svg(file)
          plot(result, cex = 0.55, print.num = F, print.pv = T, cex.pv = 0.55)
          dev.off()

        }
      )

    })



  })




}

