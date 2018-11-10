defaultspread[q_] :=  Max /@(Abs@ArrayReshape[
  Table[(q[[co, x]] - q[[co, y]])/
	(q[[co + numsymptomvariants, x]] + q[[co + numsymptomvariants, y]])
      , {co, numsymptomvariants}
      , {x, numsnpvariants - 1}
      , {y, x + 1, numsnpvariants}]
				 , {numsymptomvariants, Binomial[numsnpvariants, 2]}
				 ]);

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
    thetas,t,theta
    },

	thetas=Array[t, numsymptomvariants];

	Table[
	  sdata=data[[;;,Join[symptoms[[symptom]],snps[[snp]]]]];

	  (* one column per snpv., one row per symptomv. *)
	  f=Table[
	    Total[Boole/@Table[z==Join[symptomvariant,snpvariant],{z,sdata}]]
	    ,{symptomvariant,symptomvariants},{snpvariant,snpvariants}];
	  

	  logprob[t_] := Block[{r2 = f + t},
			       Total@Flatten@LogGamma[r2] -
			       Total@LogGamma[Total@r2] +
			       numsnpvariants*(LogGamma[Total@t] -
					       Total@LogGamma[t]) +
			       Total[logpriortheta[t]]
			 ];

	  theta=FindArgMax[{logprob[thetas], Sequence@@((# > 0)& /@thetas)}, T[{thetas, Mean@T@f}]];

	  newtheta=f+theta;
	  newA=Total@newtheta;
	  quantities = {
	    (* EV *)
	    Sequence@@T[T[newtheta]/newA],
	    (* STD *)
	    Sequence@@Sqrt[T[((newA - T@newtheta)*T@newtheta)/(newA^2*(1 + newA))]],
	   (* (* Quantiles *)
	    Sequence @@ (T@
			  Table[Quantile[
			    BetaDistribution[Sequence @@ (Reverse[newtheta[[;; , al]]])], 
			    quantiles], {al, binoutcomes + 1}]), *)
	    (* f + theta *)
	    Sequence@@newtheta,
	    (* posterior theta *)
	    Sequence@@T[Join[{theta}, Table[Null, numsnpvariants - 1, numsymptomvariants]]]
	  };

	  Join[
	    quantities,
	    T[Join[{spread[quantities]}, 
		   Table[Null, numsnpvariants - 1, numsymptomvariants]]]
		     ]

	,{symptom,numsymptoms},{snp,numsnps}]
  ];
