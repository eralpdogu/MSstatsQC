#' A function to test random forest classifiers for QC data
#' @param Test.set comma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each run.
#' @param guide.set comma-separated (.csv), metric file. It should contain a "Precursor" column and the metrics columns. It should also include "Annotations" for each run.
#' @param rf_model the model that was trained previosly by MSstatsQC-ML training process
#' @return Probability of failure predictions based on a trained model
#' @keywords machine learning, deployment
#' @export
#' @import ggplot2
#' @import FrF2
#' @importFrom dplyr recode combine
#' @import reshape2
#' @importFrom plyr ddply ldply
#' @importFrom stats median mad runif sd
#' @importFrom car boxCox yjPower
#' @importFrom h2o h2o.init as.h2o h2o.splitFrame h2o.assign h2o.grid h2o.getGrid h2o.getModel h2o.performance h2o.predict h2o.confusionMatrix
#' @examples
#' S9Site54.dataML <- DataProcess(MSstatsQC::S9Site54[, ])
#' colnames(S9Site54.dataML)[1] <- c("idfile")
#' colnames(S9Site54.dataML)[2] <- c("peptide")
#' S9Site54.dataML$peptide <- as.factor(S9Site54.dataML$peptide)
#' S9Site54.dataML$idfile <- as.numeric(S9Site54.dataML$idfile)
#' S9Site54.dataML <- within(S9Site54.dataML, rm(Annotations, missing))
#' guide.set <- dplyr::filter(S9Site54.dataML, idfile <= 20)
#' \donttest{
#' rf_model <- MSstatsQC.ML.trainR(guide.set, sim.size = 10)
#' }
#' Test.set <- dplyr::filter(S9Site54.dataML, idfile > 20)
#' \donttest{
#' MSstatsQC.ML.deployR(Test.set, guide.set, rf_model = rf_model)
#' }
#' Test.set <- S9Site54.dataML
#' \donttest{
#' MSstatsQC.ML.deployR(Test.set, guide.set, rf_model = rf_model)
#' }
#'
MSstatsQC.ML.deployR <- function(Test.set, guide.set, rf_model) {
    guide.set$peptide <- as.factor(guide.set$peptide)
    ############################################################################

    Results <- list()
    Results_annotated <- list()
    Test.set.features <- list()
    interpret.plots <- list()

    for (i in 1:nlevels(Test.set$peptide)) {
        Test.set.scale <- Test.set[Test.set$peptide == levels(Test.set$peptide)[i], c(1, 3:(ncol(Test.set)))]

        guide.set.new <- guide.set[guide.set$peptide == levels(guide.set$peptide)[i], c(3:(ncol(guide.set)))]

        for (k in 2:ncol(Test.set.scale)) {
            Test.set.scale[, k] <- (Test.set.scale[, k] - median(guide.set.new[, (k - 1)])) / mad(guide.set.new[, (k - 1)])
        }

        guide.set.new <- robust.scale(guide.set.new)

        for (k in 2:ncol(Test.set.scale)) {
            Test.set.scale[, k] <- bctrans.test((guide.set.new[, k - 1]), Test.set.scale[, k])
        }

        names(Test.set.scale) <- colnames(Test.set[, c(1, 3:(ncol(Test.set)))])

        Test.set.scale.temp <- add_features(Test.set.scale[, 2:ncol(Test.set.scale)])
        # Test.set.scale.temp <- Test.set.scale[,2:ncol(Test.set.scale)]
        Test.set.scale.temp <- Test.set.scale.temp[, order(names(Test.set.scale.temp), decreasing = TRUE)]

        Test.set.scale.h2o <- as.h2o(Test.set.scale.temp)

        Predict <- as.data.frame(h2o.predict(rf_model, Test.set.scale.h2o, type = "prob"))

        Predict <- cbind(idfile = Test.set.scale[, 1], Predict)
        Results[[i]] <- Predict[, c(1, 3)]
        Results_annotated[[i]] <- Predict$predict
        # colnames(Results)[i]<-levels(Test.set$peptide)[i]
        # colnames(Results_annotated)[i]<-levels(Test.set$peptide)[i]

        Test.set.features[[i]] <- cbind(Test.set.scale.temp, idfile = 1:length(Test.set.scale[, 1]))
        Test.set.features[[i]] <- melt(as.data.frame(Test.set.features[[i]]), id.vars = "idfile")
        g0 <- eval(substitute(
            ggplot(Test.set.features[[i]][-1, ], aes(idfile, variable)) +
                geom_tile(aes(fill = value), colour = "white") +
                labs(x = "Time", y = NULL) +
                removeGrid() +
                scale_y_discrete(expand = c(0, 0)) +
                scale_fill_gradient(
                    low = "white",
                    high = "darkorange",
                    # limits=c(-15, 100),
                    # breaks=c(0,50,100),
                    name = "Standardized and\nengineered feature values"
                ) +
                ggtitle(label = levels(Test.set$peptide)[i]) +
                theme(
                    legend.title = element_text(size = 8), legend.key.size = unit(0.5, "cm"),
                    legend.key.height = unit(0.5, "cm"), legend.justification = "bottom",
                    legend.position = "bottom", panel.background = element_blank(),
                    plot.background = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
                    axis.ticks.length = unit(0, "pt")
                ),
            list(i = i)
        ))
        interpret.plots[[i]] <- g0
    }

    FAIL <- NA
    id <- cbind(Test.set[Test.set$peptide == levels(Test.set$peptide)[which.max(table(Test.set$peptide))], 1], FAIL)
    colnames(id) <- c("idfile", "FAIL")

    for (i in 1:nlevels(Test.set$peptide)) {
        Results[[i]] <- dplyr::left_join(as.data.frame(id), Results[[i]], by = c("idfile" = "idfile"))
        Results[[i]] <- Results[[i]] %>% distinct()
        Results[[i]] <- Results[[i]][, c(1, 3)]
        colnames(Results[[i]]) <- c("idfile", "FAIL")
    }

    Results.new <- Results[[1]]
    for (i in 2:nlevels(Test.set$peptide)) {
        Results.new <- dplyr::full_join(Results.new, Results[[i]], by = "idfile")
    }

    colnames(Results.new) <- c("idfile", levels(Test.set$peptide))
    Results <- data.frame(RUN = 1:(dim(Results.new)[1]), Results.new[, -1])
    Results_melt <- melt(Results[-1, ], id.vars = "RUN")
    decision.map <- ggplot(Results_melt, aes(RUN, variable)) +
        geom_tile(aes(fill = value), colour = "white") +
        labs(x = "Time", y = NULL) +
        removeGrid() +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_gradient(
            low = "white",
            high = "red",
            limits = c(0, 1),
            breaks = c(0, 0.5, 1),
            name = "Probability\nof failure"
        ) +
        theme(
            legend.title = element_text(size = 8), legend.key.size = unit(0.5, "cm"),
            legend.key.height = unit(0.5, "cm"), legend.justification = "bottom",
            legend.position = "bottom", panel.background = element_blank(),
            plot.background = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"),
            axis.ticks.length = unit(0, "pt")
        )

    interpret.plots

    c("Drew the plots for interpretation")

    decision.map

    c("Drew the plot for final evaluation")
}
