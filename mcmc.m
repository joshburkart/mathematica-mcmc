(* ::Package:: *)

(* Mathematica Package *)

BeginPackage["MCMC`", {"Units`"}];

MCMCModelFit::usage = "MCMCModelFit[data, errors, model, paramspec, ivars, numsteps]

1. data must be given as:
{{ivar1, dvar1}, {ivar2, dvar2}, ..., {ivarN, dvarN}}
where ivar is the independent variable, dvar is the dependent variable, and N is the number of data points. Either the ivars or dvars can be vector valued; if the independent variable is a vector, then we're just dealing with a function of multiple variables, and if the dependent variable is, then we have a vector field.

2. errors must have the same length as data (= N), with each entry giving the errors in the corresponding dependent variable supplied in data. If each dvar is just a number, then so too should be each element of errors; if instead each dvar is a vector, then each element of errors should also be a vector of the same length.

3. model should evaluate to either a number or a numerical vector (depending on dvar) when all parameters and independent variables are set.

4. paramspec either gives the results of a previous MCMC run (w/ same model, data, etc.--just to add on more iterations), or lists the model parameters like so:
{{param1, ival1, spread1, domain1}, ...}
a) Each param should be symbolic.
b) ival is the initial parameter value.
c) spread is roughly how far to try to change the parameter each step in the Markov chain. In this routine we select new parameters values based on an exponential distribution of the form Exp[\[CapitalDelta]param/spread]. My Numerical Recipes book advises setting these spreads so that the average candidate acceptance is 10-40%.
d) Each domain is either Reals or a list of all possible values the parameter can take on (needs to be a uniform grid).

5. ivars gives a list of symbolic independent variables, in the same order as in data, on which model depends. If there's only one, then it need not be a list.

6. numsteps is the number of Markov chain steps to perform.";
MCMCResult::usage = "If your result object is named mcmcobj, try: mcmcobj[\"Properties\"].";

Clear["MCMC`*"];

Begin["`Private`"];

(*ChisqPlog[x_?NumericQ, degfree_?NumericQ] :=
	ChisqPlogCompiled[x, degfree];

ChisqPlogCompiled = Compile[{{x, _Real}, {degfree, _Integer}},
	-degfree/2*Log[2.] - x/2 + (degfree/2 - 1) Log[x] - Log[Gamma[degfree/2]]];*)

DiscExpNorm = Compile[{{i, _Integer}, {NN, _Integer}, {t, _Real}},
	(1 + Exp[1/t] - Exp[(i - NN)/t] - Exp[(1 - i)/t])/(Exp[1/t] - 1)];

DiscExpNormList = Compile[{{i, _Integer, 1}, {NN, _Integer, 1}, {t, _Real, 1}},
	(1 + Exp[1/t] - Exp[(i - NN)/t] - Exp[(1 - i)/t])/(Exp[1/t] - 1)];

(*DiscExpP = Compile[{{i, _Integer}, {j, _Integer}, {NN, _Integer}, {t, _Real}},
	Exp[-Abs[j - i]/t]/DiscExpNorm[i, NN, t]];*)

DiscExpPlog = Compile[{{i, _Integer}, {j, _Integer}, {NN, _Integer}, {t, _Real}},
	-Abs[j - i]/t - Log[DiscExpNorm[i, NN, t]]];

DiscExpSample = Compile[{{i, _Integer}, {NN, _Integer}, {t, _Real}, {alpha, _Real}},
		Max[If[DiscExpNorm[i, NN, t] alpha <= Exp[1/t] (1 - Exp[-i/t])/(Exp[1/t] - 1),
			Ceiling[t Log[1 + DiscExpNorm[i, NN, t] alpha (Exp[1/t] - 1) Exp[(i - 1)/t]]]
		,
			Ceiling[i - t Log[DiscExpNorm[i, NN, t] alpha (1 - Exp[1/t]) + Exp[1/t] + 1 - Exp[-(i - 1)/t]]]
		], 1]
	];

DiscExpSampleList[i_List, NN_List, t_List, alpha_List] := DiscExpSample @@@ Transpose[{i, NN, t, alpha}];

ExpSample = Compile[{state, spreads, alpha},
	state - Sign[1/2 - alpha] spreads *	(Log[2.] + Log[Min[alpha, 1 - alpha]])
];

ExpSampleList[state_List, spreads_List, alpha_List] := ExpSample @@@ Transpose[{state, spreads, alpha}];

Chisq[dpoints_List, modpoints_List, errors_List] /; Length[Dimensions[modpoints]] == 2 :=
	Total[Flatten[(modpoints - dpoints)^2 / errors^2]];
GetChisqExpr[data_List, errors_List, model_, vars_List] :=
	Module[{ipoints, dpoints, modpoints, modfunc},
		If[NumericQ[data[[1, 1]]],
			ipoints = List /@ data[[All, 1]],
			ipoints = data[[All, 1]]
		];

		If[NumericQ[data[[1, 2]]],
			dpoints = List /@ data[[All, 2]];
			modfunc = Function[Evaluate[vars], {Evaluate[model]}];
		,
			dpoints = data[[All, 2]];
			modfunc = Function[Evaluate[vars], Evaluate[model]];
		];

		modpoints = modfunc @@@ ipoints;

		Chisq[dpoints, modpoints, errors]
	];

TimeLeft[timesofar_, fractiondone_] := If[fractiondone == 0., 60 * 60 * 24. - 1., timesofar * (1. / fractiondone - 1.)];

Clear[TimeProgress];
TimeProgress[timesofar_?NumericQ, fractiondone_?NumericQ] :=
	Row[{ProgressIndicator[fractiondone],
		", Time elapsed: " <> DateString[timesofar, {"Hour24", ":", "Minute", ":", "Second"}],
		", Time left: " <>	DateString[TimeLeft[timesofar, fractiondone], {"Hour24", ":", "Minute", ":", "Second"}]}];

Sp[x__List] /; (Equal @@ Length /@ {x}) && Length[{x}] > 1 :=
		Transpose[{x}];

(*Gets y[i+1] - y[i]*)
GetDifferences[list_List] :=
		Drop[(RotateLeft[list] - list), -1];

TestMCMCMFInput[data_List, errors_List, model_, paramspec_, vars_, num_Integer] :=
(
	If[!MatchQ[paramspec, _MCMCResult],
		If[!(Length[Dimensions[paramspec]] == 2 && Dimensions[paramspec][[2]] == 4 &&
			MatchQ[paramspec[[All, 1]], {__Symbol}]),
			Message[MCMCModelFit::badinp, "bad parameter specification"];
			Return[False]
		]
	];

	If[num < 2,
		Message[MCMCModelFit::badinp, "need at least 2 steps"];
		Return[False]
	];

	If[!(MatchQ[data, {{{__?NumericQ}, {__?NumericQ}}..}] ||
		MatchQ[data, {{_?NumericQ, {__?NumericQ}}..}] ||
		MatchQ[data, {{_?NumericQ, _?NumericQ}..}] ||
		MatchQ[data, {{{__?NumericQ}, _?NumericQ}..}]),

		Message[MCMCModelFit::badinp, "data shaped inconsistently/incorrectly"];
		Return[False]
	];

	If[Length[data /. _?NumericQ -> 1 // Union] > 1,
		Message[MCMCModelFit::badinp, "data shaped inconsistently"];
		Return[False]
	];

	If[!(data[[All, 2]] /. _?NumericQ -> 1) === (errors /. _?NumericQ -> 1),
		Message[MCMCModelFit::badinp, "data shaped differently than errors"];
		Return[False]
	];

	If[!If[# == 0, 1, #]&[Length[data[[1,1]]]] == If[# == 0, 1, #]&[Length[vars]],
		Message[MCMCModelFit::badinp, "# of independent vars in data different from specified"];
		Return[False]
	];

	(*If[Head[model] === List,
		If[!Length[data[[1,2]]] == Length[model],
			Message[MCMCModelFit::badinp, "# of dependent vars in data different from model"];
			Return[False]
		]
	,
		If[!NumericQ[data[[1,2]]],
			Message[MCMCModelFit::badinp, "# of dependent vars in data different from model"];
			Return[False]
		];
	];*)

	Return[True];
);

Clear[MCMCModelFit];
Options[MCMCModelFit] = {
	"BurnFraction" -> 0.1,
	"Debug" -> False,
	"ProgressMonitor" -> Column[{Row[{"Step", "/", "MaxSteps", "  ", TimeProgress["TimeElapsed", "DoneFraction"]}],
		"CurrentParameters"(*,
		Row[{"Average acceptance: ", "AverageAcceptance"}]*)}],
	"ProgressInterval" -> 10,
	"SaveTo" -> None,
	"SaveInterval" -> 1000,
	"MakeBestFitPlot" -> False(*,
	"Compiled" -> False*)
};

MCMCModelFit::nonnumer = "Log probability given supplied model does not evaluate to a number for initial parameters; instead evaluated to: `1`\nAbort!";
MCMCModelFit::badinp = "Bad input: `1`.";
MCMCModelFit[data_List, errors_List, model_, paramspec_, invars_, num_Integer, opts : OptionsPattern[]] /;
	TestMCMCMFInput[data, errors, model, paramspec, invars, num] :=
	Module[{chisq, plog, params, spreads, state, stateval, sets, Ns, discrete,
		continuous, stateplog, candplog, cand, candval, hist, prevhist, prevnum,
		prevtime, n, i, t1, t2, burn, alpha, transplog, status = "Initializing...", resume,
		bestfitparams, z, corr, vars, plogfunc, resultlist, bestfitplot},

		plogfunc := (-#1 / 2 &);

		If[Length[invars] == 0,
			vars = {invars},
			vars = invars
		];

		Monitor[
		If[Head[paramspec] === List, (*is user attempting to resume previous mcmc run?*)
			resume = False; (*no*)

			params = paramspec[[All, 1]];
			stateval = paramspec[[All, 2]];
			spreads = paramspec[[All, 3]];
			sets = paramspec[[All, 4]];

			prevhist = {};
			prevnum = 0;
			prevtime = 0;
		,
			resume = True; (*yes*)

			params = "Parameters" /. paramspec[[1]];
			stateval = Last["ParameterRun" /. paramspec[[1]]];
			spreads = "ProposalSpreads" /. paramspec[[1]];
			sets = "ParameterDomains" /. paramspec[[1]];

			prevhist := Drop[Sp["ParameterRun", "ParametersLogPRun", "TransitionLogPRun"] /. paramspec[[1]], -1];
			prevnum = Length[prevhist];
			prevtime = "TimeSpent" /. paramspec[[1]]
		];

		n = Length[spreads]; (* # of parameters *)
		Ns = Length /@ sets; (* size of each parameter's domain; 0 for real-valued parameters *)

		status = "Evaluating chisq...";
		chisq = GetChisqExpr[data, errors, model, vars];
		plog = If[False(*OptionValue["Compiled"]*),
			status = "Compiling chisq...";
			Compile[Evaluate[params], Evaluate[plogfunc[chisq, Times @@ Dimensions[data[[All, 2]]] - n]], CompilationTarget -> "C"]
		,
			Function[Evaluate[params], Evaluate[plogfunc[chisq, Times @@ Dimensions[data[[All, 2]]] - n]]]
		];

		discrete = Flatten[Position[sets, _List]]; (* list of parameters that are discrete valued *)
		continuous = Complement[Range[n], discrete]; (* same, but instead continuous valued *)

		cand = state = stateval; (* set initial condition *)
		If[! discrete === {},
			(* ensure discrete parameters' ICs are within their domains *)
			state[[discrete]] = Nearest[sets[[#]] -> Range[Length[sets[[#]]]], stateval[[#]]][[1]] & /@ discrete;
			(* convert discrete parameters' proposal dist spreads to index spreads. doesn't really work for nonuniform domains... *)
			spreads[[discrete]] /= (Mean[GetDifferences[sets[[#]]]] & /@ discrete);
		];

		t1 = t2 = AbsoluteTime[];

		status = "Initial step...";
		If[resume,
			stateplog = Last["ParametersLogPRun" /. paramspec[[1]]];
		,
			stateplog = plog @@ stateval;
		];
		candplog = 0;

		If[!NumericQ[stateplog],
			Message[MCMCModelFit::nonnumer, stateplog];
			Return[$Failed];
		];

		(* hist is complete run history. set initial point. *)
		hist = Table[{0., 0., 0.}, {num}];

		For[i = 2, i <= num, i++,
			hist[[i-1]] = {stateval, stateplog, transplog(*Min[0, candplog - stateplog]*)};

			(* save all information about run at intervals, if desired *)
			If[Head[OptionValue["SaveTo"]] === String && Mod[i, OptionValue["SaveInterval"]] == 0,
				Put[{
						"BestFitParameters" -> Rule @@@ Sp[params, Mean[Join[prevhist, hist][[1 ;; i - 1 + prevnum, 1]]] // N],
						"ParameterErrors" -> Rule @@@ Sp[params, StandardDeviation[Join[prevhist, hist][[1 ;; i - 1 + prevnum, 1]]] // N],
						"AverageAcceptance" -> N[Mean[Exp[Join[prevhist, hist][[1 ;; i - 1 + prevnum, 3]]]]],
						"TimeSpent" -> (t2 - t1) Second + prevtime,
						"Parameters" -> params,
						"ProposalSpreads" -> spreads,
						"ParameterDomains" -> sets,
						"BurnFraction" -> OptionValue["BurnFraction"],
						"BurnEnd" -> burn,
						"ParameterRun" -> Join[prevhist, hist][[1 ;; i - 1 + prevnum, 1]],
						"ParametersLogPRun" -> Join[prevhist, hist][[1 ;; i - 1 + prevnum, 2]],
						"TransitionLogPRun" -> Join[prevhist, hist][[1 ;; i - 1 + prevnum, 3]],
						"BurnFraction" -> OptionValue["BurnFraction"]
					}
				,
					OptionValue["SaveTo"]
				]
			];

			(* update status indicator *)
			If[Mod[i, OptionValue["ProgressInterval"]] == 0,
				status = OptionValue["ProgressMonitor"] /. \
					{
						"CurrentParameters" -> Rule @@@ Sp[params, stateval],
						"TimeElapsed" -> t2 - t1,
						"DoneFraction" -> (i - 1)/(num),
						"Step" -> i,
						"MaxSteps" -> num,
						"AverageAcceptance" -> If[And[! FreeQ[OptionValue["ProgressMonitor"], "AverageAcceptance"], i > 2],
							Chop[N[Mean[Exp[hist[[2 ;; i-1, 3]]]]]], Null
						],
						"Plot" -> If[OptionValue["MakeBestFitPlot"] && Length[vars] == 1
							&& And[! FreeQ[OptionValue["ProgressMonitor"], "Plot"], i > 2],
							Show[
								Plot[Evaluate[model /. (Rule @@@ Sp[params, stateval]) /. vars[[1]] -> z], {z, Min[data[[All, 1]]], Max[data[[All, 1]]]}],
								If[Length[data[[1, 2]]] == 0,
									ListPlot[data, Joined->False],
									ListPlot[Table[Sp[data[[All, 1]], data[[All, 2, i]]], {i, 1, Length[model]}], Joined->False]
								],
								Frame -> True,
								Axes -> False
							]
						,
							"Number of ind. variables > 1."
						]
					}
			];

			alpha = RandomReal[{0, 1}, n]; (* random variables with which to generate candidate point *)
			If[! continuous === {},
				cand[[continuous]] = ExpSampleList[state[[continuous]], spreads[[continuous]], alpha[[continuous]]];
			];
			If[! discrete === {},
				cand[[discrete]] = DiscExpSampleList[state[[discrete]], Ns[[discrete]], spreads[[discrete]], alpha[[discrete]]]
			];

			candval = cand;
			If[! discrete === {},
				candval[[discrete]] = sets[[#, cand[[#]]]] & /@ discrete;
			];

			If[OptionValue["Debug"], Print["cand ",cand," ","candval",candval]];

			candplog = plog @@ candval;

			transplog = Min[0., candplog - stateplog +
				If[! discrete === {}, (* discrete proposal dist is not symmetric (continuous is). take this into account. *)
					Total[DiscExpPlog @@@ Sp[state[[discrete]], cand[[discrete]], Ns[[discrete]], spreads[[discrete]]]] -
						Total[DiscExpPlog @@@ Sp[cand[[discrete]], state[[discrete]], Ns[[discrete]], spreads[[discrete]]]]
				,
					0
				]
			];

			alpha = RandomReal[{0, 1}]; (* random variable with which to determine whether to accept candidate *)

			If[OptionValue["Debug"], Print["transplog ",transplog," vs. random plog ",Log[alpha]]];

			If[Log[alpha] < transplog,
				 If[OptionValue["Debug"], Print["cand accepted; had plog ",transplog]];

				 state = cand;
				 stateval = candval;
				 stateplog = candplog
				 ,
				 If[OptionValue["Debug"], Print["cand rejected; had plog ",transplog]];
			];
			t2 = AbsoluteTime[];
		];

		hist[[num]] = {stateval, stateplog, Min[0, candplog - stateplog]};

		burn = Ceiling[Min[(num + prevnum)/2, Max[1000, (num + prevnum)*OptionValue["BurnFraction"]]]];

		bestfitparams = Rule @@@ Sp[params, Mean[Join[prevhist, hist][[burn ;; num + prevnum, 1]]] // N];

		status = "Computing correlation matrix...";
		corr = Correlation[Join[prevhist, hist][[All, 1]]];

		status = "Done!";
		, status];

		resultlist = {
			"BestFitParameters" -> bestfitparams,
			"ParameterErrors" -> Rule @@@ Sp[params, StandardDeviation[Join[prevhist, hist][[burn ;; num + prevnum, 1]]] // N],
			"AverageAcceptance" -> N[Mean[Exp[Join[prevhist, hist][[burn ;; num + prevnum, 3]]]]],
			"TimeSpent" -> (t2 - t1) Second + prevtime,
			"NumSteps" -> num + prevnum,
			"Parameters" -> params,
			"ProposalSpreads" -> spreads,
			"ParameterDomains" -> sets,
			"BurnFraction" -> OptionValue["BurnFraction"],
			"BurnEnd" -> burn,
			"BestFitReducedChisq" -> (chisq /. bestfitparams) / (Length[data] - n),
			"CorrelationMatrix" -> MatrixForm[corr],
			"ParameterRun" -> Join[prevhist, hist][[All, 1]],
			"ParametersLogPRun" -> Join[prevhist, hist][[All, 2]],
			"TransitionLogPRun" -> Join[prevhist, hist][[All, 3]]
		};
		If[OptionValue["MakeBestFitPlot"],
			bestfitplot = If[Length[vars] == 1,
				Show[
					Plot[Evaluate[model /. bestfitparams /. vars[[1]] -> z], {z, Min[data[[All, 1]]], Max[data[[All, 1]]]}, PlotRange -> All],
					If[Length[data[[1, 2]]] == 0,
						ListPlot[data, Joined->False, PlotRange -> All],
						ListPlot[Table[Sp[data[[All, 1]], data[[All, 2, i]]], {i, 1, Length[model]}], Joined->False, PlotRange -> All]
					],
					Frame -> True,
					Axes -> False
				]
			,
				"Number of ind. variables > 1."
			];
			resultlist = Append[resultlist, "BestFitPlot" -> bestfitplot];
		];
		MCMCResult[resultlist]
	];

Clear[MCMCResult];
Format[MCMCResult[list_List]] := "MCMCResult"["BestFitParameters" /. list, "\[LeftSkeleton]" <> ToString[Length["ParameterRun" /. list]] <> "\[RightSkeleton]"];

(this:MCMCResult[list_List])["ParameterRunPlots", opts___] :=
	Table[
		ListPlot[Transpose[("ParameterRun" /. list)][[i]],
			AxesLabel -> {"Step", ToString[this["Parameters"][[i]]]},
			FrameLabel -> {"Step", ToString[this["Parameters"][[i]]]},
			opts
		],
	{i, Length[this["Parameters"]]}];

(this:MCMCResult[list_List])["ParameterHistograms", opts___] :=
	Table[If[this["NumSteps"] > 1*^6,
		SmoothHistogram[#,
			Filling -> Axis,
			Axes -> {True, False},
			Ticks -> {Automatic, None},
			AxesLabel -> {ToString[this["Parameters"][[i]]], None},
			opts]&
	,
		Histogram[#,
			Ticks -> {Automatic, None},
			Axes -> {True, False},
			AxesLabel -> {ToString[this["Parameters"][[i]]], None},
			opts]&][
		Transpose[("ParameterRun" /. list)[[this["BurnEnd"] ;; this["NumSteps"]]]][[i]]
	], {i, Length[this["Parameters"]]}];

(this:MCMCResult[list_List])["Properties"] := Join[list[[All, 1]], {"ParameterRunPlots", "ParameterHistograms"}];

(this:MCMCResult[list_List])[str_String] := str /. list;
(this:MCMCResult[list_List])[{str__String}] := Rule @@@ Sp[{str}, this /@ {str}];

End[];

EndPackage[];
