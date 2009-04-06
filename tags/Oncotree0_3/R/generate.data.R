"generate.data" <-
function(N, otree, with.errors=TRUE,
          edge.weights=if (with.errors) "estimated" else "observed"){
	distr <- distribution.oncotree(otree, with.probs=TRUE, with.errors=with.errors,
	      edge.weights=edge.weights)
	epos <- otree$eps["epos"]
	eneg <- otree$eps["eneg"]
	ran.idx <- sample(1:nrow(distr), size=N, prob=distr$Prob, replace=TRUE)
	ran.data <- distr[ran.idx, 2:otree$nmut]
	rownames(ran.data) <- 1:N
	ran.data
}

