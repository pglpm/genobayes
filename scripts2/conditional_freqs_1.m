condfreqstatistics[data_,
		       symptoms_,
		       symptomvariants_,
		       snps_,
		       snpvariants_,
		       namesymptoms_,
		       namesymptomvariants_,
		       namesnps_,
		       namesnpvariants_,
		       savedir_,
		       filename_,
		       logpriortheta_,
		       spread_,
		       writeflag_:True] :=
  Block[{
    numsymptoms = Length[symptoms],
    numsymptomvariants = Length[symptomvariants],
    numsnps = Length[snps],
    numsnpvariants = Length[snpvariants],
    namestatistics = Flatten@Table[ii<>jj,{ii,{"EV_", "SD_", "post.theta_", "opt.theta_","max.spread_"}},{jj,namesymptomvariants}],
    thetas
    },

	thetas=Unique["t",numsymptomvariants];

	Table[
	  sdata=data[[;;,{symptoms[[symptom]],snps[[snp]]}]];

	  f=Table[
	    Total[Boole/@Table[z==Join[symptomvariant,snpvariant],{z,sdata}]]
	    ,{snpvariant,snpvariants},{symptomvariant,symptomvariants}];

	  logprob[t_] := Block[{r2 = f + t},
			       Total@Flatten@LogGamma[r2] -
			       Total@LogGamma[Total@r2] +
			       noutcomes*(LogGamma[Total@t] - Total@LogGamma[t]) +
			       Total[logpriortheta[Log@t]]
			 ];

	  theta=thetas/.FindMaximum[{logprob[thetas], Sequence@@((# > 0)& /@thetas)}, T[{thetas, Mean[f]}]][[2]];



	,{symptom,numsymptoms},{snp,numsnps}]
  ]
