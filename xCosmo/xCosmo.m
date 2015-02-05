(* ::Package:: *)

(* ::Title:: *)
(*xCosmo*)


(* ::Text:: *)
(*:Author:*)
(*Joseph Elliston*)


(* ::Text:: *)
(*:Summary:*)
(*This package simplifies the use of xAct in performing various common calculations in cosmology. *)


(* ::Text:: *)
(*:Context: *)
(*xAct`xCosmo`*)


(* ::Text:: *)
(*:Copyright:*)
(*Copyright (C) University of Sussex 2013-2014*)


(* ::Text:: *)
(*:History:*)
(*V1.0 - 27.08.2014: First release*)


(* ::Section::Closed:: *)
(*Preliminaries*)


(* ::Subsection::Closed:: *)
(*GPL*)


(* ::Text:: *)
(*Copyright (C) 2013-2014 University of Sussex*)
(**)
(*This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.*)
(**)
(*This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.*)
(**)
(*You should have received a copy of the GNU General Public License along with xAct; if not, write to the Free Software Foundation, Inc., 59 Temple Place-Suite 330, Boston, MA 02111-1307, USA. *)


(* ::Subsection::Closed:: *)
(*Version numbers and package dependencies*)


(* ::Text:: *)
(*Package version number:*)


xAct`xCosmo`$Version={"1.0",Date[][[{1,2,3}]]};


(* ::Text:: *)
(*Expected version number for xTensor:*)


xAct`xCosmo`$xTensorVersionExpected={"1.0.5",{2013,1,30}};


(* ::Text:: *)
(*In case the package is loaded multiple times, wipe the memory of all package symbols apart from the version numbers defined above.*)


With[
	{xAct`xCosmo`Private`TensorDecompositionSymbols=
		DeleteCases[
			Join[Names["xAct`xCosmo`*"],Names["xAct`xCosmo`Private`*"]],
			"$Version"|"xAct`xCosmo`$Version"|"$xTensorVersionExpected"|"xAct`xCosmo`$xTensorVersionExpected"
		]
	},
	Unprotect/@xAct`xCosmo`Private`TensorDecompositionSymbols;
	Clear/@xAct`xCosmo`Private`TensorDecompositionSymbols;
]


(* ::Text:: *)
(*Set this package to be the last package to load*)


If[Unevaluated[xAct`xCore`Private`$LastPackage]===xAct`xCore`Private`$LastPackage,
	xAct`xCore`Private`$LastPackage="xAct`xCosmo`"
];


(* ::Text:: *)
(*Begin the package and load dependencies*)


BeginPackage["xAct`xCosmo`",{"xAct`xCore`","xAct`xPerm`","xAct`xTensor`","xAct`xPert`","xAct`SplitExpression`","xAct`ByParts`"}];


(* ::Text:: *)
(*Check version of xTensor:*)


If[
	Not@OrderedQ@Map[Last,{$xTensorVersionExpected,xAct`xTensor`$Version}],
	Message[General::versions,"xTensor",xAct`xTensor`$Version,$xTensorVersionExpected];
	Abort[]
];


(* ::Subsection::Closed:: *)
(*Output message after loading the package*)


Print[xAct`xCore`Private`bars];
Print["Package xAct`xCosmo version ",$Version[[1]],", ",$Version[[2]]];
Print["Copyright (C) 2013-2014, University of Sussex, under the General Public License."];
Print["Written by Joseph Elliston."];


(* ::Text:: *)
(*We specify the context xAct`xCosmo` to avoid overriding the Disclaimer of xCore, xPerm and xTensor. However we need to turn off the message General:shdw temporarily:*)


Off[General::shdw];
xAct`xCosmo`Disclaimer[]:=Print["These are points 11 and 12 of the General Public License:\n\nBECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM `AS IS\.b4 WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n\nIN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES."];
On[General::shdw];


(* ::Text:: *)
(*If xAct`xCosmo` is the last package read, then print the short GPL disclaimer:*)


If[xAct`xCore`Private`$LastPackage==="xAct`xCosmo`",
Unset[xAct`xCore`Private`$LastPackage];
Print[xAct`xCore`Private`bars];
Print["These packages come with ABSOLUTELY NO WARRANTY; for details type Disclaimer[]. This is free software, and you are welcome to redistribute it under certain conditions. See the General Public License for details."];
Print[xAct`xCore`Private`bars]];


(* ::Subsection::Closed:: *)
(*Reduce and configure output*)


(* ::Text:: *)
(*Some of this is already done by the SplitExpression package.*)


(* ::Text:: *)
(*The default for System`$Assumptions is True. Reset this to {} so that we can easily append new assumptions.*)


$Assumptions = {};


(* ::Subsection::Closed:: *)
(*Detailed package description*)


(* ::Text:: *)
(*Please see the file 'Pedagogical guide' for a detailed package description*)


(* ::Subsection::Closed:: *)
(*Acknowledgements*)


(* ::Text:: *)
(*My thanks go to David Seery and the University of Sussex for allowing me to devote so much time into writing this package. *)
(*I am also indebted to Guido Pettinari for patiently sharing some of his vast knowledge of Mathematica.*)
(*--Joe Elliston*)


(* ::Section::Closed:: *)
(*Public context*)


(* ::Subsection::Closed:: *)
(*Usage messages*)


(* ::Subsubsection::Closed:: *)
(*Control variables*)


Global`SingleScalarField::usage="If True then xCosmo will presume that there is a single scalar field \[CurlyPhi] and a scalar potental V is also defined as a scalar function. If SingleScalarField=False, then \[CurlyPhi] is defined as a rank-1 tensor with a field space index. A field space manifold is also defined, along with a field metrc G. In the multi-field case the potential V is a scalar tensor.";
Global`decomposition::usage = "Determines which decomposition of the spacetime metric shall be used. Options are ADM, FlatFLRW and UniformDensityFLRW. Additional spatial metrics can easily be added by editing the xCosmo.m file.";


(* ::Subsubsection::Closed:: *)
(*Functions*)


Unprotect[FieldEquation];
FieldEquation::usage = "FieldEquation[\[Delta]tensor] computes the field equations from the Lagrangian \[ScriptCapitalL] which comes from considering the variation of \[ScriptCapitalL] with respect to the tensor variation \[Delta]tensor. The function internally expands \[ScriptCapitalL] to first order in perturbations and then uses VarD to pull out the result. If \[Delta]tensor has indices, such as \[Delta]g[LI[1],\[Mu],\[Nu]] then the indices of the resulting equation will be complementary to these indices (covariant in this case).";

Unprotect[ExpandAllOrdersQuantities];
ExpandAllOrdersQuantities::usage = "ExpandAllOrdersQuantities[order][expr] acts on expr to expand any all-orders quantities up to the order 'order' and then returns the result at that prescribed order, not including lower order pieces.";

Unprotect[PerturbativeExpansionRules];
PerturbativeExpansionRules::usage = "PerturbativeExpansionRules[order] are the rules defined internally by xCosmo to define how quantities are to be expanded perturbatively for a given decompositions of the spacetime metric g.";

Unprotect[hsubs];
hSubs::usage ="hSubs[expr] substitutes for the 3-metric h and its determinant in expr.";


(* ::Subsubsection::Closed:: *)
(*Spacetime quantities*)


Unprotect[MSpaceTime];
MSpaceTime::usage=
"The spacetime manifold, with metric g and uses lowercase Greek indices \[Mu],\[Nu],\[Gamma],\[Eta],\[Kappa],\[Lambda],\[Rho],\[Sigma],\[Theta],\[Tau],\[Upsilon],\[Chi],\[Omega]. Other indices can be added using the function AddIndices[TangentMSpaceTime,{new indices}].";

Unprotect[TangentMSpaceTime];
TangentMSpaceTime::usage="Tangent Vbundle to the MSpaceTime manifold. The greek spacetime indices are associated to this VBundle.";

Unprotect[g];
g::usage="g[\[Mu],\[Nu]] is the 4-dimensional spacetime metric.";

Unprotect[CDSpaceTime];
CDSpaceTime::usage=
"Covariant derivative associated to the spacetime metric g.";

Unprotect[\[CurlyEpsilon]];
\[CurlyEpsilon]::usage="Constant parameter to denote perturbative orders.";

Unprotect[\[Delta]g];
\[Delta]g::usage="\[Delta]g[LI[order],\[Mu],\[Nu]] is the perturbation to the perturbed 4-dimensional spacetime metric. In xCosmo this is only used internally by the FieldEquation function.";

Unprotect[Global`MassPlanck];
Global`MassPlanck::usage = "The constant Planck mass. The global $Assumptions is configured to allow Mathematica to know that MassPlanck>0.";



(* ::Subsubsection::Closed:: *)
(*Lagrangian*)


\[ScriptCapitalL]::usage = "Reserved symbol in xCosmo for use as the scalar Lagrangian. This is read in by the functions FieldEquation.";

Unprotect[\[CurlyPhi]];
\[CurlyPhi]::usage = "The full Scalar field \[CurlyPhi][] in the single field mode (or \[CurlyPhi][A] in the multi-field mode). \[CurlyPhi] = background + perturbations, and therefore is spatially dependent.";

Unprotect[\[Phi]];
\[Phi]::usage = "Background scalar field, written as \[Phi][] in the single field mode or \[Phi][A] in the multi-field mode. This only depends on time.";

Unprotect[\[Delta]\[Phi]];
\[Delta]\[Phi]::usage = "Perturbation to \[CurlyPhi][] or \[CurlyPhi][A] in single or multi-field mode respectively. Use as \[Delta]\[Phi][LI[order]] in the single field mode or \[Delta]\[Phi][LI[order],A] in the multi-field mode.";

Unprotect[V];
V::usage = "The scalar field potential, assumed to only be a function of spatially-dependent \[CurlyPhi] and not its time derivatives. In the single field case this is a scalar function defined via DefScalarFunction. In the multi-field case this is a scalar Tensor on the MFieldSpace manifold, with an implicit dependence on \[CurlyPhi]. The background potential in the single field case is written as V[\[Phi]] and in the multifield case as Vbar[].";


(* ::Subsubsection::Closed:: *)
(*Lap and InvLap*)


Unprotect[UseLaplacian];
UseLaplacian::usage = "If UseLaplacian=True then the functions MakeLap and ExpandLap are activated, otherwise if UseLaplacian=False then these functions do nothing.";

Unprotect[VBundleOfLaplacian];
VBundleOfLaplacian::usage = "Defines the Vbundle that indices have to be a member of in order that they can be reformulated by the Laplacian functions.";

Unprotect[Lap];
Lap::usage = "Spatial Laplacian, Lap[X] = \[Delta][i,j]PD[-i]@PD[-j]@X. Lap is defined to be distributive over terms added together. To determine which terms are inert under the Laplacian function, we state that any term that is not inert must be dependent on MSpace. Any other terms are removed from the Lap argument. For example, the scalar \[CurlyPhi][] is spatially dependent and so to signify this we set it to be depend on MSpace. Lap is an inert head, meaning that index validation can see through this head.";

Unprotect[InvLap];
InvLap::usage = "Inverse spatial Laplacian. InvLap is defined to be distributive over terms added together. Inert terms are also removed, as explained in the Lap function usage message. InvLap is also an an inert head.";

Unprotect[MakeLap];
MakeLap::usage = "MakeLap[arg] finds instances of the spatial Laplacian in 'arg' and writes them in terms of the Lap function.";

Unprotect[ExpandLap];
ExpandLap::usage = "Expands out any factors of Lap[X] into index notation as Lap[X] \[Rule] \[Delta][i,j]PD[-i]@PD[-j]@X.";


(* ::Subsection::Closed:: *)
(*Spacetime manifolds and metrics*)


(* ::Text:: *)
(*Define the Planck Mass as a constant, and add its positive definite attribute to the list of internal assumptions.*)


DefConstantSymbol[Global`MassPlanck,PrintAs->"\!\(\*SubscriptBox[\(M\), \(Pl\)]\)"];
AppendTo[$Assumptions,Global`MassPlanck>0];


(* ::Text:: *)
(*Define the spacetime manifold and indices. The available indices for spacetime are all of the lowercase Greek letters listed below.*)


DefManifold[MSpaceTime,4,{\[Mu],\[Nu],\[Gamma],\[Eta],\[Kappa],\[Lambda],\[Rho],\[Sigma],\[Theta],\[Tau],\[Upsilon],\[Chi],\[Omega]}];


(* ::Text:: *)
(*Define the spacetime metric 'g' and associated covariant derivative CDSpaceTime. We also define \[Delta]g as the metric perturbations, and \[CurlyEpsilon] as the parameter that counts perturbative orders. The metric perturbations \[Delta]g are only required for the FieldEquation function.*)


DefMetric[-1,g[-\[Mu],-\[Nu]],CDSpaceTime,{";","\[Del]"}];
DefMetricPerturbation[g,\[Delta]g,\[CurlyEpsilon]];


(* ::Text:: *)
(*Now tell Mathematica that $PerturbationParameter == \[CurlyEpsilon] is a constant and should not be differentiated :*)


SetAttributes[{$PerturbationParameter, \[CurlyEpsilon]}, Constant];


(* ::Text:: *)
(*Make curvature quantities output in more pleasant ways:*)


PrintAs[RiemannCDSpaceTime]^="R";
PrintAs[RicciCDSpaceTime]^="R";
PrintAs[RicciScalarCDSpaceTime]^="R";
PrintAs[EinsteinCDSpaceTime]^="G";
PrintAs[RiemannCDSpace]^="R";
PrintAs[RicciCDSpace]^="R";
PrintAs[RicciScalarCDSpace]^="\!\(\*SuperscriptBox[\(R\), \((3)\)]\)";
PrintAs[EinsteinCDSpace]^="G";


(* ::Subsection::Closed:: *)
(*Define single scalar field*)


(* ::Text:: *)
(*In this case Global`SingleScalarField = True.*)
(*\[CurlyPhi] is the full scalar field, \[Phi] is the background field and \[Delta]\[Phi] is the perturbation.*)
(*V is a scalar function that we shall use as the scalar potential.*)
(**)
(*A general tensor X is expanded as X -> X[LI[0]] + \[Delta]X[LI[1]] + \[Delta]X[LI[2]] + ... where we do not include factors of 1/2, 1/6 etc.*)
(*The only exception to this expansion rule is for the scalar field \[Phi][] (the fields \[Phi][A] in the multi-field case). In this case there is no obligation to perturbations beyond linear order, since we can define \[Delta]\[Phi][LI[1]] or \[Delta]\[Phi][LI[1],A] as the total perturbation to all orders. Doing so yields simpler analysis. Expanding into different perturbative orders can be easily done if, for example, one may wish to define \[Delta]\[Phi][LI[1]] or \[Delta]\[Phi][LI[1],A] as the Gaussian perturbation and higher order perturbations then denote the non-Gaussianity. *)


If[Global`SingleScalarField,
	DefTensor[\[CurlyPhi][],{MSpaceTime}];
	DefTensor[\[Phi][],{}];
	DefTensorPerturbation[\[Delta]\[Phi][LI[order]],\[CurlyPhi][],{}];
	DefScalarFunction[V];

	(*Define the perturbation scheme*)
	\[Delta]\[Phi]/:\[Delta]\[Phi][LI[n_]]:=0/;n>1;
	Unprotect[Perturbation];
	Perturbation/:HoldPattern[Perturbation[V[B_],n_Integer]]/;(n>1):=Perturbation[\[Delta]\[Phi][LI[1]]D[V[B],B],n-1];
	Perturbation/:HoldPattern[Perturbation@V[B_]]:=\[Delta]\[Phi][LI[1]]D[V[B],B];
	Perturbation/:HoldPattern[Perturbation[A__@V[B_],n_Integer]]/;(n>1):=Perturbation[\[Delta]\[Phi][LI[1]]D[A@V[B],B],n-1];
	Perturbation/:HoldPattern[Perturbation@A__@V[B_]]:=\[Delta]\[Phi][LI[1]]D[A@V[B],B];
	Protect[Perturbation];
];


(* ::Subsection::Closed:: *)
(*Define multiple scalar fields and a field-metric G*)


(* ::Text:: *)
(*If Global`SingleScalarField = False, define the manifold MFieldSpace and associated metric G and covariant derivative CDFieldSpace. Nominally, this is a 2-dimensional manifold, but the dimensionality is irrelevant unless we want to specify the field metric G. Indices in the tangent space of MFieldSpace are uppercase Roman. Other letters we have left out because they are either special Mathematica characters or they are useful elsewhere, e.g. we use 'H' as the Hubble rate.*)
(**)
(*xCosmo is currently configured for a flat field metric G and so we remove any terms involving derivatives of this metric. *)
(*The perturbations are also themselves covariantly constant.*)


If[!Global`SingleScalarField,
Unprotect[MFieldSpace];
MFieldSpace::usage=
"(Only defined if SingleScalarField = False) The manifold on which the multiple scalar fields \[CurlyPhi][A] are defined, known as the field-space manifold. MFieldSpace has metric G which is currently presumed to be Euclidean and uses uppercase Roman indices A,B,F,J,L,M,P,Q,R,S. Other indices can be added using the function AddIndices[TangentMFieldSpace,{new indices}].";
Unprotect[TangentMFieldSpace];
TangentMFieldSpace::usage="(Only defined if SingleScalarField = False) Tangent bundle to the MFieldSpace manifold.";
DefManifold[MFieldSpace,2,{A,B,F,J,L,M,P,Q,R,S}];

Unprotect[G];
G::usage="(Only defined if SingleScalarField = False) Metric of the field space on which the scalar fields \[CurlyPhi][A] are defined. Currently assumed to be Euclidean.";
Unprotect[CDFieldSpace];
CDFieldSpace::usage="(Only defined if SingleScalarField = False) Covariant derivative associated to the field-space metric G.";
DefMetric[-1,G[-A, -B],CDFieldSpace,{";", "\[Del]"},FlatMetric->True];

DefTensor[\[CurlyPhi][A],MFieldSpace];
DefTensor[\[Phi][A],MFieldSpace];
DefTensorPerturbation[\[Delta]\[Phi][LI[order],A],\[CurlyPhi][A],MFieldSpace];

(*Define the perturbation scheme*)
\[Delta]\[Phi]/:\[Delta]\[Phi][LI[n_],_]:=0/;n>1;
G/:CD_[_]@G[__]:=0/;CovDQ[CD];
\[Delta]\[Phi]/:CDFieldSpace[_]@\[Delta]\[Phi][__]:=0;
Unprotect[Perturbation];
Perturbation/:HoldPattern[Perturbation[V[],n_Integer]]/;(n>1):=Module[{A},Perturbation[\[Delta]\[Phi][LI[1],A]CDFieldSpace[-A]@V[],n-1]];
Perturbation/:Perturbation@V[]:=Module[{A},\[Delta]\[Phi][LI[1],A]CDFieldSpace[-A]@V[]];
Perturbation/:HoldPattern[Perturbation[stuff__@V[],n_Integer]]/;(n>1):=Module[{A},Perturbation[\[Delta]\[Phi][LI[1],A]CDFieldSpace[-A]@stuff@V[],n-1]];
Perturbation/:HoldPattern[Perturbation@stuff__@V[]]:=Module[{A},\[Delta]\[Phi][LI[1],A]CDFieldSpace[-A]@stuff@V[]];
Perturbation/:Perturbation[G[__],n_Integer]:=0;
Perturbation/:Perturbation@G[__]:=0;
Protect[Perturbation];

DefTensor[V[],{}];
Unprotect[Vbar];
Vbar::usage = "(Only defined if SingleScalarField = False) In the multifield mode, xCosmo uses V[] for the full potential and Vbar[] for the background potential. This distinction is not needed in the single field mode because V is then a scalar function and this can be written as a scalar function of \[CurlyPhi] or \[Phi] to differentiate the two cases.";
DefTensor[Vbar[],{},PrintAs->"\!\(\*OverscriptBox[\(V\), \(_\)]\)"];
(*The following line of code ensures that when VarD is varying V then it expands V out in terms of \[CurlyPhi]
We include a similar result for Vbar, but this should never actually be used.*)
V/:ImplicitTensorDepQ[V,\[CurlyPhi]]=True;
Vbar/:ImplicitTensorDepQ[Vbar,\[Phi]]=True;

(*Define how to take parametric derivatives of V and Vbar*)
(*ParamD/:ParamD[a_,b___]@A_:=Module[{B},ParamD[b]@CDFieldSpace[-B]@A ParamD[a]@\[CurlyPhi][B]]
	/;(Position[A,V[],Infinity]=!={})/;(Last@InverseComposition@A===V[]);
ParamD/:ParamD[a_,b___]@A_:=Module[{B},ParamD[b]@CDFieldSpace[-B]@A ParamD[a]@\[Phi][B]]
	/;(Position[A,Vbar[],Infinity]=!={})/;(Last@InverseComposition@A===Vbar[]);
*)

PrintAs[RiemannCDFieldSpace]^="R";
PrintAs[RicciCDFieldSpace]^="R";
PrintAs[RicciScalarCDFieldSpace]^="R";
PrintAs[EinsteinCDFieldSpace]^="G";
];


(* ::Subsection::Closed:: *)
(*Decompositions*)


(* ::Subsubsection::Closed:: *)
(*Preliminary discussion*)


(* ::Text:: *)
(*There are many possible ways to decompose teh spacetime metric and each choice may in general require a different splitting of the 4 spacetime dimensions, e.g. 3+1, 2+1+1, 1+1+1+1. The control variable decomposition tells xCosmo which of the decompositions below to use. Each decomposition has a separate section below enclosed in an If function of the form If[decomposition=MyDecomposition, <<definitions>>]. Within these If functions we do the following:*)


(* ::Item:: *)
(*Define the manifolds, metrics and parameters required by the splitting that shall be used.*)


(* ::Item:: *)
(*Ensure that any existing tensors have the correct paramter dependence*)


(* ::Item:: *)
(*Prescrive if the Laplacian functions will be used*)


(* ::Item:: *)
(*Define the tensors that will exist after the index splitting*)


(* ::Item:: *)
(*State any assumptions that need to be known about these new tensors*)


(* ::Item:: *)
(*Define the splitting rules*)


(* ::Item:: *)
(*State how any all-orders quantities will be expanded out into individual perturbative orders*)


(* ::Item:: *)
(*Define rules to automatically account for the special ways that derivatives can act on tensors after the splitting, e.g. imposing the spatial independence of some tensors.*)


(* ::Text:: *)
(*We also define any usage rules at the same time, so that the user is only provided with usage rules that are relevant to their choice of decomposition.*)


(* ::Text:: *)
(*Note on the identity matrix: Some metrics are naturally written in terms of the identity matrix. This can provide a challenge because the object that raises and lowers indices is in most cases not the Euclidean metric. Errors could then occur if the indices on the identity tensor were moved up/down. Brief experimentation with defining a second 'frozen' metric did not seem fruitful because index contractions appear to not be automatic (or I was doing it wrong...). Instead we have a robust solution by defining separate tensors tensors, \[Delta]up and \[Delta]down for use in the different situations. Should the indices of these tensors be raised or lowered then we will automatically keep track of where they should be and their proper placement can always be reinstated by using the SeparateMetric function.*)


(* ::Subsubsection::Closed:: *)
(*ADM*)


(* ::Text:: *)
(*This is the classic 3+1 split into Lapse and Shift where we do not manipulate the spatial metric further. We define a spatial manifold and a Time parameter. The available spatial indices the lowercase Roman letters from i to z. Define also the spatial metric 'h' with associated covariant derivative CDSpace. *)
(**)
(*The Lapse and Shiftare both all-orders quantities and the shift is zero at background order. Individual perturbative orders of the Lapse and scalar Shift are respectively denoted by the tensors \[Alpha] and \[Beta]. The individual orders of the pure-vector shift are denoted by the tensor \[Beta]vec, which we make divergenceless. The background order \[Alpha] and \[Beta] and \[Beta]vec are known from the outset.*)
(**)
(*This is also where we set terms to be zero under the action of a PD derivative. This is used later in the Lap and InvLap functions to determine which terms are inert.*)


If[Global`decomposition==="ADM",

	(******************Define manifolds, parameters, metrics and tensors******************)

	Unprotect[MSpace];
	MSpace::usage="The 3-space manifold, with metric h and uses lowercase Roman indices running alphabetically from i through to z. Other indices can be added using the function AddIndices[TangentMSpace,{new indices}].";
	Unprotect[TangentMSpace];
	TangentMSpace::usage="Tangent bundle to the MSpace manifold.";
	DefManifold[MSpace,3,IndexRange[i,z]];
	
	Unprotect[h];
	h::usage="h[i,j] is the 3-dimensional spatial metric.";
	Unprotect[CDSpace];
	CDSpace::usage="Covariant derivative associated to the spacial 3-metric h.";
	DefMetric[1,h[-i,-j],CDSpace,{"|","\[ScriptCapitalD]"}];
	
	Time::usage="Time is a parameter.";
	DefParameter[Time,PrintAs->"t"];
	

	(******************Add any parameter dependencies to existing tensors as appropriate******************)

	AddTensorDependencies[{h,g,\[Delta]g,\[CurlyPhi],\[Phi],\[Delta]\[Phi]},{Time}];
	If[!Global`SingleScalarField,
		AddTensorDependencies[{V,Vbar},{Time}];
	];

	(******************Define usage of the Laplacian functions******************)
	UseLaplacian=True;
	VBundleOfLaplacian=TangentMSpace;

	(******************Define the splitting and the splitting quantities******************)

	DefSplitting[ADM,MSpaceTime->{Time,MSpace}];

	Unprotect[Lapse];
	Lapse::usage="The lapse scalar to all perturbative orders, printed as N.";
	DefTensor[Lapse[],{Time},PrintAs->"N"];

	Unprotect[Shift];
	Shift::usage="The shift vector to all perturbative orders, printed as N with an index.";
	DefTensor[Shift[i],{Time,MSpace},PrintAs->"N"];

	Unprotect[\[Alpha]];
	\[Alpha]::usage="The lapse scalar perturbation.";
	DefTensor[\[Alpha][LI[order]],{Time}];

	Unprotect[\[Beta]];
	\[Beta]::usage="The shift scalar perturbation.";
	DefTensor[\[Beta][LI[order]],{Time}];

	Unprotect[\[Beta]vec];
	\[Beta]vec::usage="The shift pure-vector perturbation, printed as \!\(\*OverscriptBox[\(\[Beta]\), \(~\)]\).";
	DefTensor[\[Beta]vec[LI[order],i],{Time,MSpace},PrintAs->"\!\(\*OverscriptBox[\(\[Beta]\), \(~\)]\)"];

	(******************Splitting rules******************)

	SplittingRules[ADM,{
	g[-Time,-Time]->-Lapse[]^2+Shift[-i]Shift[i],
	g[-Time,-i]->Shift[-i],
	g[-i,-Time]->Shift[-i],
	g[-i,-j]->h[-i,-j],

	g[Time,Time]->-Lapse[]^-2,
	g[Time,i]->Lapse[]^-2 Shift[i],
	g[i,Time]->Lapse[]^-2 Shift[i],
	g[i,j]->h[i,j]-Lapse[]^-2 Shift[i]Shift[j],

	Detg[]->-Lapse[]^2 Deth[]
	}];

	hSubs=#&;

	(******************Join any assumptions for these new tensors to the global $Assumptions list******************)
	$Assumptions=Join[$Assumptions,{Lapse[]>0}];

	(******************Rules for expansion of all-orders quantities******************)

	PerturbativeExpansionRules[order_]:={
		Lapse[]:>1 + Sum[\[CurlyEpsilon]^ii \[Alpha][LI[ii]],{ii,1,order}],
		Shift[a_]:>Sum[\[CurlyEpsilon]^ii(Module[{j},h[a,j]PD[-j]@\[Beta][LI[ii]]]+\[Beta]vec[LI[ii],a]),{ii,1,order}]
	};

	(******************Define rules for how any tensors behave under any particular derivatives******************)

	Unprotect[PD];
	(*The traceless nature of the vector perturbation*)
	PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Beta]vec[_,_],Infinity]=!={})/;MatchQ[Last@InverseComposition@A,\[Beta]vec[_,-a]];
	(*Spatially independent quantities*)
	If[Global`SingleScalarField,
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Phi][],Infinity]=!={})/;(Last@InverseComposition@A===\[Phi][]);
		,
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Phi][_],Infinity]=!={})/;MatchQ[Last@InverseComposition@A,\[Phi][_]];
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},Vbar[],Infinity]=!={})/;(Last@InverseComposition@A===Vbar[]);
	];	
	Protect[PD];
];


(* ::Subsubsection::Closed:: *)
(*FlatFLRW*)


(* ::Text:: *)
(*Following the ADM split, the spatial metric is then prescribed to be of the flat form.*)


If[Global`decomposition==="FlatFLRW",

	(******************Define manifolds, parameterS and metrics******************)

	Unprotect[MSpace];
	MSpace::usage="The 3-space manifold, with metric h and uses lowercase Roman indices running alphabetically from i through to z. Other indices can be added using the function AddIndices[TangentMSpace,{new indices}].";
	Unprotect[TangentMSpace];
	TangentMSpace::usage="Tangent bundle to the MSpace manifold.";
	DefManifold[MSpace,3,IndexRange[i,z]];
	
	Unprotect[h];
	h::usage="h[i,j] is the 3-dimensional spatial metric.";
	Unprotect[CDSpace];
	CDSpace::usage="Covariant derivative associated to the spacial 3-metric h.";
	DefMetric[1,h[-i,-j],CDSpace,{"|","\[ScriptCapitalD]"}];
	
	Time::usage="Time is a parameter.";
	DefParameter[Time,PrintAs->"t"];	

	(******************Add any parameter dependencies to existing tensors as appropriate******************)

	AddTensorDependencies[{h,g,\[Delta]g,\[CurlyPhi],\[Phi],\[Delta]\[Phi]},{Time}];
	If[!Global`SingleScalarField,
		AddTensorDependencies[{V,Vbar},{Time}];
	];

	(******************Define usage of the Laplacian functions******************)
	UseLaplacian=True;
	VBundleOfLaplacian=TangentMSpace;

	(******************Define the splitting and the splitting quantities******************)

	DefSplitting[FlatFLRW,MSpaceTime->{Time,MSpace}];

	Unprotect[Lapse];
	Lapse::usage="The lapse scalar to all perturbative orders, printed as N.";
	DefTensor[Lapse[],{Time},PrintAs->"N"];

	Unprotect[Shift];
	Shift::usage="The shift vector to all perturbative orders, printed as N with an index.";
	DefTensor[Shift[i],{Time,MSpace},PrintAs->"N"];

	Unprotect[\[Alpha]];
	\[Alpha]::usage="The lapse scalar perturbation.";
	DefTensor[\[Alpha][LI[order]],{Time}];

	Unprotect[\[Beta]];
	\[Beta]::usage="The shift scalar perturbation.";
	DefTensor[\[Beta][LI[order]],{Time}];

	Unprotect[\[Beta]vec];
	\[Beta]vec::usage="The shift pure-vector perturbation, printed as \!\(\*OverscriptBox[\(\[Beta]\), \(~\)]\).";
	DefTensor[\[Beta]vec[LI[order],i],{Time,MSpace},PrintAs->"\!\(\*OverscriptBox[\(\[Beta]\), \(~\)]\)"];

	Unprotect[\[Delta]up];
	\[Delta]up::usage = "\[Delta]up[i,j] is the tensor equivalent to the contravariant 3-dimensional metric on a Eucldean space.";
	DefTensor[\[Delta]up[i,j],MSpace,Symmetric[{i,j}],PrintAs->"\[Delta]"];
	
	Unprotect[\[Delta]down];
	\[Delta]down::usage = "\[Delta]down[-i,-j] is the tensor equivalent to the covariant 3-dimensional metric on a Eucldean space.";
	DefTensor[\[Delta]down[-i,-j],MSpace,Symmetric[{-i,-j}],PrintAs->"\[Delta]"];

	Unprotect[a];
	a::usage = "The FLRW scale factor, only a function of the parameter Time.";
	DefTensor[a[],Time];

	Unprotect[H];
	H::usage = "The FLRW Hubble rate, only a function of the parameter Time.";
	DefTensor[H[],Time];

	(******************Join any assumptions for these new tensors to the global $Assumptions list******************)
	$Assumptions=Join[$Assumptions,{a[]>0,Lapse[]>0}];

	(******************Configure the basic rules for \[Delta]up and \[Delta]down******************)

	\[Delta]up/:A_[__][\[Delta]up[__]]:=0/;(CovDQ[A]||A===ParamD||A===LieD);
	\[Delta]down/:A_[__][\[Delta]down[__]]:=0/;(CovDQ[A]||A===ParamD||A===LieD);

	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[i,-l]/;(k===j)/;!MemberQ[ABIndexQ/@{j,k},False];
	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[i,-k]/;(l===j)/;!MemberQ[ABIndexQ/@{j,l},False];
	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[j,-l]/;(k===i)/;!MemberQ[ABIndexQ/@{i,k},False];
	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[j,-k]/;(l===i)/;!MemberQ[ABIndexQ/@{i,l},False];

	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[i,-l]/;(k===j)/;!MemberQ[ABIndexQ/@{j,k},False];
	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[i,-k]/;(l===j)/;!MemberQ[ABIndexQ/@{j,l},False];
	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[j,-l]/;(k===i)/;!MemberQ[ABIndexQ/@{i,k},False];
	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[j,-k]/;(l===i)/;!MemberQ[ABIndexQ/@{i,l},False];

	(******************Splitting rules******************)

	SplittingRules[FlatFLRW,{
	g[-Time,-Time]->-Lapse[]^2+Shift[-i]Shift[i],
	g[-Time,-i]->Shift[-i],
	g[-i,-Time]->Shift[-i],
	g[-i,-j]->h[-i,-j],

	g[Time,Time]->-Lapse[]^-2,
	g[Time,i]->Lapse[]^-2 Shift[i],
	g[i,Time]->Lapse[]^-2 Shift[i],
	g[i,j]->h[i,j]-Lapse[]^-2 Shift[i]Shift[j],

	Detg[]->-Lapse[]^2 a[]^6
	}];

	hSubs=#//.{
		h[i_,j_]/;(ABIndexQ@i&&ABIndexQ@j)/;(UpIndexQ@i&&UpIndexQ@j):>a[]^-2 \[Delta]up[i,j],
		h[i_,j_]/;(ABIndexQ@i&&ABIndexQ@j)/;(DownIndexQ@i&&DownIndexQ@j):>a[]^2 \[Delta]down[i,j],
		Deth[]->a[]^6
	}&;

	(******************Rules for expansion of all-orders quantities******************)

	PerturbativeExpansionRules[order_]:={
		Lapse[]:>1 + Sum[\[CurlyEpsilon]^ii \[Alpha][LI[ii]],{ii,1,order}],
		Shift[a_]:>Sum[\[CurlyEpsilon]^ii(Module[{j},h[a,j]PD[-j]@\[Beta][LI[ii]]]+\[Beta]vec[LI[ii],a]),{ii,1,order}]
	};

	(******************Define rules for how any tensors behave under any particular derivatives******************)

	(*Time derivatives of the scale factor are written in terms of the Hubble rate*)
	a/:ParamD[Time]@a[]:=a[]H[];
	a/:ParamD[Time,Time]@a[]:=a[](ParamD[Time]@H[]+H[]^2);

	Unprotect[PD];
	(*The traceless nature of the vector perturbation*)
	PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Beta]vec[_,_],Infinity]=!={})/;MatchQ[Last@InverseComposition@A,\[Beta]vec[_,-a]];
	(*Spatially independent quantities*)
	PD/:PD[b_]@A_:=0/;(VBundleOfIndex@b===TangentMSpace)/;(Position[{A},a[],Infinity]=!={})/;(Last@InverseComposition@A===a[]);
	PD/:PD[b_]@A_:=0/;(VBundleOfIndex@b===TangentMSpace)/;(Position[{A},H[],Infinity]=!={})/;(Last@InverseComposition@A===H[]);
	If[Global`SingleScalarField,
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Phi][],Infinity]=!={})/;(Last@InverseComposition@A===\[Phi][]);
		,
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Phi][_],Infinity]=!={})/;MatchQ[Last@InverseComposition@A,\[Phi][_]];
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},Vbar[],Infinity]=!={})/;(Last@InverseComposition@A===Vbar[]);
	];	
	Protect[PD];
];


(* ::Subsubsection::Closed:: *)
(*UniformDensityFLRW*)


(* ::Text:: *)
(*Following the ADM split, the spatial metric is then prescribed to be of the Uniform density form, which includes the \[Zeta] scalar. This degree of freedom may then be removed by prescribing that the density is constant at each perturbative order. This constraint is, however, not automatically included.*)


If[Global`decomposition==="UniformDensityFLRW",

	(******************Define manifolds, parameterS and metrics******************)

	Unprotect[MSpace];
	MSpace::usage="The 3-space manifold, with metric h and uses lowercase Roman indices running alphabetically from i through to z. Other indices can be added using the function AddIndices[TangentMSpace,{new indices}].";
	Unprotect[TangentMSpace];
	TangentMSpace::usage="Tangent bundle to the MSpace manifold.";
	DefManifold[MSpace,3,IndexRange[i,z]];
	
	Unprotect[h];
	h::usage="h[i,j] is the 3-dimensional spatial metric.";
	Unprotect[CDSpace];
	CDSpace::usage="Covariant derivative associated to the spacial 3-metric h.";
	DefMetric[1,h[-i,-j],CDSpace,{"|","\[ScriptCapitalD]"}];
	
	Time::usage="Time is a parameter.";
	DefParameter[Time,PrintAs->"t"];	

	(******************Add any parameter dependencies to existing tensors as appropriate******************)

	AddTensorDependencies[{h,g,\[Delta]g,\[CurlyPhi],\[Phi],\[Delta]\[Phi]},{Time}];
	If[!Global`SingleScalarField,
		AddTensorDependencies[{V,Vbar},{Time}];
	];

	(******************Define usage of the Laplacian functions******************)
	UseLaplacian=True;
	VBundleOfLaplacian=TangentMSpace;

	(******************Define the splitting and the splitting quantities******************)

	DefSplitting[UniformDensityFLRW,MSpaceTime->{Time,MSpace}];

	Unprotect[Lapse];
	Lapse::usage="The lapse scalar to all perturbative orders, printed as N.";
	DefTensor[Lapse[],{Time},PrintAs->"N"];

	Unprotect[Shift];
	Shift::usage="The shift vector to all perturbative orders, printed as N with an index.";
	DefTensor[Shift[i],{Time,MSpace},PrintAs->"N"];

	Unprotect[\[Alpha]];
	\[Alpha]::usage="The lapse scalar perturbation.";
	DefTensor[\[Alpha][LI[order]],{Time}];

	Unprotect[\[Beta]];
	\[Beta]::usage="The shift scalar perturbation.";
	DefTensor[\[Beta][LI[order]],{Time}];

	Unprotect[\[Beta]vec];
	\[Beta]vec::usage="The shift pure-vector perturbation, printed as \!\(\*OverscriptBox[\(\[Beta]\), \(~\)]\).";
	DefTensor[\[Beta]vec[LI[order],i],{Time,MSpace},PrintAs->"\!\(\*OverscriptBox[\(\[Beta]\), \(~\)]\)"];

	Unprotect[\[Delta]up];
	\[Delta]up::usage = "\[Delta]up[i,j] is the tensor equivalent to the contravariant 3-dimensional metric on a Eucldean space.";
	DefTensor[\[Delta]up[i,j],MSpace,Symmetric[{i,j}],PrintAs->"\[Delta]"];
	
	Unprotect[\[Delta]down];
	\[Delta]down::usage = "\[Delta]down[-i,-j] is the tensor equivalent to the covariant 3-dimensional metric on a Eucldean space.";
	DefTensor[\[Delta]down[-i,-j],MSpace,Symmetric[{-i,-j}],PrintAs->"\[Delta]"];

	Unprotect[a];
	a::usage = "The FLRW scale factor, only a function of the parameter Time.";
	DefTensor[a[],Time];

	Unprotect[H];
	H::usage = "The FLRW Hubble rate, only a function of the parameter Time.";
	DefTensor[H[],Time];

	Unprotect[\[Zeta]];
	\[Zeta]::usage = "The curvature perturbation on uniform density hypersurfaces at a give perturbative order.";
	DefTensor[\[Zeta][LI[order]],{Time,MSpace}];

	Unprotect[zeta];
	zeta::usage = "The all-orders curvature perturbation on uniform density hypersurfaces.";
	DefTensor[zeta[],{Time,MSpace},PrintAs->"\[Zeta]"];


	(******************Join any assumptions for these new tensors to the global $Assumptions list******************)
	$Assumptions=Join[$Assumptions,{a[]>0,Lapse[]>0}];

	(******************Configure the basic rules for \[Delta]up and \[Delta]down******************)

	\[Delta]up/:A_[__][\[Delta]up[__]]:=0/;(CovDQ[A]||A===ParamD||A===LieD);
	\[Delta]down/:A_[__][\[Delta]down[__]]:=0/;(CovDQ[A]||A===ParamD||A===LieD);

	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[i,-l]/;(k===j)/;!MemberQ[ABIndexQ/@{j,k},False];
	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[i,-k]/;(l===j)/;!MemberQ[ABIndexQ/@{j,l},False];
	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[j,-l]/;(k===i)/;!MemberQ[ABIndexQ/@{i,k},False];
	\[Delta]up/:\[Delta]up[i_,j_]*\[Delta]down[-k_,-l_]:=h[j,-k]/;(l===i)/;!MemberQ[ABIndexQ/@{i,l},False];

	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[i,-l]/;(k===j)/;!MemberQ[ABIndexQ/@{j,k},False];
	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[i,-k]/;(l===j)/;!MemberQ[ABIndexQ/@{j,l},False];
	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[j,-l]/;(k===i)/;!MemberQ[ABIndexQ/@{i,k},False];
	\[Delta]up/:\[Delta]down[-k_,-l_]*\[Delta]up[i_,j_]:=h[j,-k]/;(l===i)/;!MemberQ[ABIndexQ/@{i,l},False];

	(******************Splitting rules******************)

	SplittingRules[UniformDensityFLRW,{
	g[-Time,-Time]->-Lapse[]^2+Shift[-i]Shift[i],
	g[-Time,-i]->Shift[-i],
	g[-i,-Time]->Shift[-i],
	g[-i,-j]->h[-i,-j],

	g[Time,Time]->-Lapse[]^-2,
	g[Time,i]->Lapse[]^-2 Shift[i],
	g[i,Time]->Lapse[]^-2 Shift[i],
	g[i,j]->h[i,j]-Lapse[]^-2 Shift[i]Shift[j],

	Detg[]->-Lapse[]^2 a[]^6 Exp[6 \[CurlyEpsilon] zeta[]]
	}];

	hSubs=#//.{
		h[i_,j_]/;(ABIndexQ@i&&ABIndexQ@j)/;(DownIndexQ@i&&DownIndexQ@j):>a[]^2 Exp[2 \[CurlyEpsilon] zeta[]] \[Delta][i,j],
		h[i_,j_]/;(ABIndexQ@i&&ABIndexQ@j)/;(UpIndexQ@i&&UpIndexQ@j):>a[]^-2 Exp[-2 \[CurlyEpsilon] zeta[]] \[Delta][i,j],
		Deth[]->a[]^6 Exp[6 \[CurlyEpsilon] zeta[]]
	}&;

	(******************Rules for expansion of all-orders quantities******************)

	PerturbativeExpansionRules[order_]:={
		Lapse[]:>1 + Sum[\[CurlyEpsilon]^ii \[Alpha][LI[ii]],{ii,1,order}],
		Shift[a_]:>Sum[\[CurlyEpsilon]^ii(Module[{j},h[a,j]PD[-j]@\[Beta][LI[ii]]]+\[Beta]vec[LI[ii],a]),{ii,1,order}],
		zeta[]:>Sum[\[CurlyEpsilon]^ii \[Zeta][LI[ii]],{ii,1,order}]
	};

	(******************Define rules for how any tensors behave under any particular derivatives******************)

	(*Time derivatives of the scale factor are written in terms of the Hubble rate*)
	a/:ParamD[Time]@a[]:=a[]H[];
	a/:ParamD[Time,Time]@a[]:=a[](ParamD[Time]@H[]+H[]^2);

	Unprotect[PD];
	(*The traceless nature of the vector perturbation*)
	PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Beta]vec[_,_],Infinity]=!={})/;MatchQ[Last@InverseComposition@A,\[Beta]vec[_,-a]];
	(*Spatially independent quantities*)
	PD/:PD[b_]@A_:=0/;(VBundleOfIndex@b===TangentMSpace)/;(Position[{A},a[],Infinity]=!={})/;(Last@InverseComposition@A===a[]);
	PD/:PD[b_]@A_:=0/;(VBundleOfIndex@b===TangentMSpace)/;(Position[{A},H[],Infinity]=!={})/;(Last@InverseComposition@A===H[]);
	If[Global`SingleScalarField,
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Phi][],Infinity]=!={})/;(Last@InverseComposition@A===\[Phi][]);
		,
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},\[Phi][_],Infinity]=!={})/;MatchQ[Last@InverseComposition@A,\[Phi][_]];
		PD/:PD[a_]@A_:=0/;(VBundleOfIndex@a===TangentMSpace)/;(Position[{A},Vbar[],Infinity]=!={})/;(Last@InverseComposition@A===Vbar[]);
	];	
	Protect[PD];
];


(* ::Section::Closed:: *)
(*Private Context*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Debugging*)


(* ::Text:: *)
(*Define the variable InfoLevel to control the output of useful debugging messages. *)


(* ::Item:: *)
(*InfoLevel=0 --- All output messages are minimised.*)


(* ::Item:: *)
(*InfoLevel=1 --- Higher-level output messages are permitted.*)


(* ::Item:: *)
(*InfoLevel=2 --- All output messages are permitted.*)


InfoLevel=0;


(* ::Text:: *)
(*Also define the function InfoPrint[stuff,level] that prints out 'stuff' only if InfoLevel is sufficient.*)


SetNumberOfArguments[InfoPrint,2];
InfoPrint[expr_,level_]:=If[level<=InfoLevel,Print[expr];];


(* ::Subsection::Closed:: *)
(*Laplacian and inverse Laplacian operators*)


(* ::Subsubsection::Closed:: *)
(*Basic rules*)


(* ::Text:: *)
(*We shall now define the Laplacian operator 'Lap' and its inverse 'InvLap'.*)
(*These are defined using InertHead to ensure that index manipulations act through these functions.*)


DefInertHead[Lap];
DefInertHead[InvLap];
LapQ=(MatchQ[#,_Lap])&;
InvLapQ=(MatchQ[#,_InvLap])&;


(* ::Text:: *)
(*Make the output format for these functions look nice:*)


PrintAs[Lap]^="\!\(\*SuperscriptBox[\(\[PartialD]\), \(2\)]\)";
PrintAs[InvLap]^="\!\(\*SuperscriptBox[\(\[PartialD]\), \(-2\)]\)";


(* ::Text:: *)
(*Make Lap and InvLap know that they are inverses of one another:*)


Lap[InvLap[aa_]] := aa;
InvLap[Lap[aa_]] := aa;


(* ::Text:: *)
(*Number rules:*)


Lap[_?NumberQ] := 0;
InvLap[_?NumberQ] := 0;


(* ::Text:: *)
(*Make Lap and InvLap distributive over terms added together (Plus heads):*)


Lap[expr_] := Lap/@Expand[expr]/;(Head@Expand@expr===Plus);
InvLap[expr_] := InvLap/@Expand[expr]/;(Head@Expand@expr===Plus);


(* ::Text:: *)
(*Pull PD and ParamD derivatives outside Lap and InvLap functions*)


Lap/:Lap[CD_[A_][B__]]:=CD[A]@Lap[B]/;(PD===CD||ParamD===CD);
InvLap/:InvLap[CD_[A_][B__]]:=CD[A]@InvLap[B]/;(PD===CD||ParamD===CD);


(* ::Subsubsection::Closed:: *)
(*Factoring out inert terms*)


(* ::Text:: *)
(*Pull inert terms out of the arguments of 'Lap'\.1d and 'InvLap'. Inert terms are any terms that are zero when differentiated by PD with a new spatial index.*)


LapInertTerms[expr_]:=
Module[{termlist,nonInertPositions,InertPositions},
	termlist=If[Head@expr===Times,List@@expr,{expr}];
	nonInertPositions=DeleteDuplicates@Flatten[First/@Position[termlist,_Lap|_InvLap|_?((xTensorQ@Head@#&&PD[-NewIndexIn@VBundleOfLaplacian]@#=!=0)&)]];
	InertPositions=Partition[DeleteCases[Array[#&,Length@termlist],Alternatives@@nonInertPositions],1];
	If[InertPositions==={},
		1,
		Times@@Extract[termlist,InertPositions]
	]
];

Lap[expr__?(LapInertTerms@#=!=1&)]:=
Module[{tmp},
	tmp=LapInertTerms@expr;
	tmp*Lap[expr/tmp]
];

InvLap[expr__?(LapInertTerms@#=!=1&)]:=
Module[{tmp},
	tmp=LapInertTerms@expr;
	tmp*InvLap[expr/tmp]
];


(* ::Subsubsection::Closed:: *)
(*'MakeLap' and 'ExpandLap' functions*)


(* ::Text:: *)
(*Define 'MakeLap' that takes an expressions and isolates Lap functions. 'ExpandLap' does the opposite.*)
(*These are defined relative to the decompostiion being used. Some decompositions, e.g. 1+1+1+1 will not have any terms of Laplacian form and in this case MakeLap and ExpandLap do nothing.*)


If[UseLaplacian,
	MakeLap=(Expand[#]//.{
		\[Delta]up[i_,j_]*terms___*PD[-k_]@PD[-l_]@A__:>terms*Lap[A]
			/;!MemberQ[VBundleOfIndex@#===VBundleOfLaplacian&/@{i,j,k,l},False]
			/;(Sort[{i,j}]===Sort[{k,l}]),
		\[Delta]up[i_,j_]*terms___*OP__@PD[-k_]@PD[-l_]@A__:>terms*OP@Lap[A]
			/;!MemberQ[VBundleOfIndex@#===VBundleOfLaplacian&/@{i,j,k,l},False]
			/;(Sort[{i,j}]===Sort[{k,l}])
	})&;
	ExpandLap=(NoScalar[
		#//.HoldPattern[Lap[X__]]:>Scalar[Module[{ii,jj},
			ii=NewIndexIn@VBundleOfLaplacian;
			jj=NewIndexIn@VBundleOfLaplacian;
			\[Delta]up[ii,jj]*PD[-jj]@PD[-ii]@X]
			]
	]//ReplaceDummies//ScreenDollarIndices)&;
	,
	MakeLap=#&;
	ExpandLap=#&;
];


(* ::Subsection::Closed:: *)
(*Field equations from the Lagrangian*)


(* ::Subsubsection::Closed:: *)
(*Varying the Lagrangian*)


(* ::Text:: *)
(*When comparing perturbative orders, xAct does not automatically understand identities involving delta of LI indices. Add these:*)


Unprotect[delta];
delta /: delta[-LI[a_], LI[b_]] := 1 /; (a == b);
delta /: delta[-LI[a_], LI[b_]] := 0 /; (a != b);
Protect[delta];


(* ::Text:: *)
(*We also teach VarD how to cope with Parametric derivatives and variations of the scalar potential V in the multifield case.*)


VarD/:VarD[A_[D___],PD][ParamD[E_,F__]@B_,C__]:=VarD[A[D],PD][ParamD[F]@B,-ParamD[E]@C];
VarD/:VarD[A_[D___],PD][ParamD[E_]@B_,C__]:=VarD[A[D],PD][B,-ParamD[E]@C];
If[!Global`SingleScalarField,
	VarD/:VarD[\[CurlyPhi][A_],_][V[],rest___]:=rest*CDFieldSpace[-A]@V[]
];


(* ::Subsubsection::Closed:: *)
(*FieldEquation function*)


(* ::Text:: *)
(*The energy momentum tensors are computed using the function FieldEquation[\[ScriptCapitalL],arg] where arg is the quantity that the Lagrangian \[ScriptCapitalL] is being varied with respect to.*)


SetNumberOfArguments[FieldEquation,1];
FieldEquation[arg_]:=
Module[{\[Delta]\[ScriptCapitalL],variation},
	\[Delta]\[ScriptCapitalL]=ExpandPerturbation@Perturbation@(Sqrt[-Detg[]]\[ScriptCapitalL])//ContractMetric//Simplification//ReplaceDummies//ScreenDollarIndices;
	variation=-2/Sqrt[-Detg[]] VarD[arg,CDSpaceTime][\[Delta]\[ScriptCapitalL]];
	If[!Global`SingleScalarField,
		variation=variation/.VarD[A_,CDSpaceTime][B_,C_]:>VarD[A,CDFieldSpace][B,C];
	];
	variation//RicciToEinstein//ContractMetric//Simplification//Expand//ScreenDollarIndices
]


(* ::Subsection::Closed:: *)
(*Perturbed actions*)


(* ::Subsubsection::Closed:: *)
(*ExpandAllOrdersQuantities*)


(* ::Input:: *)
(*?ExpandAllOrdersQuantities*)


ExpandAllOrdersQuantities[order_][expr_]:=
Module[{output},
	output=FixedPoint[ExpandAllOrdersQuantitiesInternal[order],expr];
	SeriesCoefficient[output,{\[CurlyEpsilon],0,order}]
];

ExpandAllOrdersQuantitiesInternal[order_][expr_]:=
Module[{newexpr},
	newexpr=expr//.PerturbativeExpansionRules[order];
	If[Global`SingleScalarField,
		newexpr=newexpr//.V[\[CurlyPhi][]]:>V[\[Phi][]]+Sum[\[CurlyEpsilon]^ii Nest[\[Delta]\[Phi][LI[1]]D[#,\[Phi][]]&,V[\[Phi][]],ii],{ii,1,order}];
		newexpr=newexpr//.\[CurlyPhi][]->\[Phi][]+\[CurlyEpsilon] \[Delta]\[Phi][LI[1]];
		,
		newexpr=newexpr//.V[]:>Vbar[]+Sum[\[CurlyEpsilon]^ii Nest[Module[{A},\[Delta]\[Phi][LI[1],A]CDFieldSpace[-A]@#]&,Vbar[],ii],{ii,1,order}];
		newexpr=newexpr//.\[CurlyPhi][A_]:>\[Phi][A]+\[CurlyEpsilon] \[Delta]\[Phi][LI[1],A]
	];
	newexpr=newexpr//hSubs;
	newexpr=Normal@Series[newexpr,{\[CurlyEpsilon],0,order}];
	newexpr//Simplification//NoScalar//MakeLap//ContractMetric//Expand
];


(* ::Section::Closed:: *)
(*Finish Up*)


(* ::Subsection::Closed:: *)
(*End Private Context*)


End[]


(* ::Subsection::Closed:: *)
(*Protect variables*)


(* ::Subsubsection::Closed:: *)
(*Functions*)


Protect[FieldEquation];
Protect[ExpandAllOrdersQuantities];
Protect[PerturbativeExpansionRules];
Protect[hsubs];


(* ::Subsubsection::Closed:: *)
(*Spacetime quantities*)


Protect[MSpaceTime];
Protect[TangentMSpaceTime];
Protect[g];
Protect[CDSpaceTime];
Protect[\[CurlyEpsilon]];
Protect[\[Delta]g];
Protect[Global`MassPlanck];


(* ::Subsubsection::Closed:: *)
(*Lagrangian*)


Protect[\[CurlyPhi]];
Protect[\[Phi]];
Protect[\[Delta]\[Phi]];
Protect[V];


(* ::Subsubsection::Closed:: *)
(*Lap and InvLap*)


Protect[UseLaplacian];
Protect[VBundleOfLaplacian];
Protect[Lap];
Protect[InvLap];
Protect[MakeLap];
Protect[ExpandLap];
Protect[hsubs];


(* ::Subsubsection::Closed:: *)
(*Scalar field*)


If[!Global`SingleScalarField,
Protect[MFieldSpace];
Protect[TangentMFieldSpace];
Protect[G];
Protect[CDFieldSpace];
Protect[Vbar];
];


(* ::Subsubsection::Closed:: *)
(*Decomposition-dependent*)


(* ::Text:: *)
(*ADM*)


If[Global`decomposition==="ADM",
Protect[MSpace];
Protect[TangentMSpace];
Protect[h];
Protect[CDSpaceTime];
Protect[CDSpace];
Protect[Lapse];
Protect[Shift];
Protect[\[Alpha]];
Protect[\[Beta]];
Protect[\[Beta]vec];
];


(* ::Text:: *)
(*FlatFLRW*)


If[Global`decomposition==="FlatFLRW",
Protect[MSpace];
Protect[TangentMSpace];
Protect[h];
Protect[CDSpaceTime];
Protect[CDSpace];
Protect[Lapse];
Protect[Shift];
Protect[\[Alpha]];
Protect[\[Beta]];
Protect[\[Beta]vec];
Protect[\[Delta]up];
Protect[\[Delta]down];
Protect[a];
Protect[H];
];


(* ::Text:: *)
(*UniformDensityFLRW*)


If[Global`decomposition==="UniformDensityFLRW",
Protect[MSpace];
Protect[TangentMSpace];
Protect[h];
Protect[CDSpaceTime];
Protect[CDSpace];
Protect[Lapse];
Protect[Shift];
Protect[\[Alpha]];
Protect[\[Beta]];
Protect[\[Beta]vec];
Protect[\[Delta]up];
Protect[\[Delta]down];
Protect[a];
Protect[H];
Protect[\[Zeta]];
Protect[zeta];
];


(* ::Subsection::Closed:: *)
(*End package*)


EndPackage[]
