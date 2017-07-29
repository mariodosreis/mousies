          seed =  -1

       seqfile = set.6s.80662 
      Imapfile = Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

 speciesdelimitation = 0 * fixed species tree
         speciestree = 0  * speciestree pSlider ExpandRatio ShrinkRatio

    speciesmodelprior = 1         * 0: uniform labeled histories; 1:uniform rooted trees; 2:user probs

  species&tree = 6    O  L  M  X  R  B
                      1  1  1  1  1  1
                     (O,((L,M),(X,(R,B))));
                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 20000  * number of data sets in seqfile

     cleandata = 1    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 2 1000   # gamma(a, b) for theta
      tauprior = 2 250 1  # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

      finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 8000
      sampfreq = 2
       nsample = 20000
