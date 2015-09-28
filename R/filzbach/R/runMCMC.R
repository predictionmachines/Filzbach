runMCMC <-
function (burnin, eststeps, loglikelihood, samplesize, paramdefs, 
    thinning = 100) 
{
    .External("filzbachRunMCMC", burnin, eststeps, loglikelihood, 
        samplesize, paramdefs, thinning, new.env())
}
