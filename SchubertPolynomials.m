BeginPackage["SchubertPolynomials`"];

Memoization::usage="Memoization is an option to SchubertPolynomial, DoubleSchubertPolynomial, GrothendieckPolynomial, DoubleGrothendieckPolynomial, SchubertExpansion, and SchubertPolynomialQ that can be set to True or False to enable or disable the use of memoization in the corresponding computation. The defualt value is True.";
Reduced::usage="Reduced is an option to Mitosis that can be set to True or False to enable or disable the inclusion of nonreduced pipe dreams in the mitosis algorithm. The default value is True.";
CrossColor::usage="CrossColor is an option to DrawPipeDream, that can be set to any color. The defualt value is Purple.";
ElbowColor::usage="ElbowColor is an option to DrawPipeDream, that can be set to any color. The defualt value is Purple.";
Labels::usage="Labels is an option to DrawPipeDream and DrawBumplessPipeDream that can be set to True or False. The defualt value is True.";
PermutationListForm::usage="PermutationListForm[w] converts the permutation w from the Cycles[] form into a list of numbers.";
PipeDreamQ::usage="PipeDreamQ[P] checks if P is a matrix of zeros and ones that is a pipe dream when zeros are viewed as elbows and ones as crossings. ";
PipeDreamToPermutation::usage = "PipeDreamToPermutation[P] returns the permutation represented by the pipe dream P.";
PipeDreams::usage = "PipeDreams[w] returns all pipe dreams (reduced and nonreduced) of the permutation w.";
DrawPipeDream::usage="DrawPipeDream[P] gives a visual representation of the pipe dream defined by the matrix P.";
Inversions::usage="Inversions[w] returns the inversion set of the permutation w, the set of all (i,j) with i<j such that w(i)>w(j).";
NumberInversions::usage = "NumberInversions[w] returns the number of inversions of the permutation w.";
ReducedPipeDreamQ::usage = "ReducedPipeDreamQ[P] returns True if P is a reduced pipe dream and False otherwise .";
ReducedPipeDreams::usage = "ReducedPipeDreams[w] returns the set of all reduced pipe dreams of the permutation w.";
GrothendieckPolynomial::usage = "GrothendieckPolynomial[w] returns the Grothendieck polynomial of w in variables x[i]. Optional flag Memoization can be set to True or False to enable automatic storing of output.";
BetaGrothendieckPolynomial::usage = "BetaGrothendieckPolynomial[w] returns the beta Grothendieck polynomial of w in variables x[i] and b.";
DoubleGrothendieckPolynomial::usage = "DoubleGrothendieckPolynomial[w] returns the double Grothendieck polynomial of w in variables x[i] and y[j]. Optional flag Memoization can be set to True or False to enable automatic storing of output.";
BetaDoubleGrothendieckPolynomial::usage="BetaDoubleGrothendieckPolynomial[w] returns the beta double Grothendieck polynomial of w in variables x[i], y[j], and b.";
SchubertPolynomial::usage = "SchubertPolynomial[w] returns the Schubert polynomial of w in variables x[i]. Optional flag Memoization can be set to True or False to enable automatic storing of output.";
Descents::usage = "Descents[w] returns the descent set of the permutation w, the set of all i such that w(i)>w(i+1).";
Ascents::usage = "Ascents[w] returns the ascent set of the permutation w, the set of all i such that w(i)<w(i+1).";
Transposition::usage = "Transposition[i,j] returns the transposition swapping i and j.";
TransitionRule::usage = "TransitionRule[w] returns {{r,v},trans(w)} where \!\(\*SubscriptBox[\(S\), \(w\)]\)=x[r] \!\(\*SubscriptBox[\(S\), \(v\)]\)+\!\(\*UnderscriptBox[\(\[Sum]\), \(trans \((w)\)\)]\)\!\(\*SubscriptBox[\(S\), \(\[Sigma]\)]\).";
DoubleSchubertPolynomial::usage = "DoubleSchubertPolynomial[w] returns the double Schubert polynomial of w in variables x[i] and y[j]. Optional flag Memoization can be set to True or False to enable automatic storing of output.";
Letter::usage ="Letter[i] returns the adjacent transposition (i i+1) (in \!\(\*SubscriptBox[\(S\), \(i + 1\)]\)).";
ReducedWords::usage ="ReducedWords[w] returns the set of all the reduced words {\!\(\*SubscriptBox[\(i\), \(1\)]\),...,\!\(\*SubscriptBox[\(i\), \(k\)]\)} such that for \!\(\*SubscriptBox[\(s\), \(j\)]\) the transposition (j j+1), w(k)=\!\(\*SubscriptBox[\(s\), SubscriptBox[\(i\), \(k\)]]\)(\!\(\*SubscriptBox[\(s\), SubscriptBox[\(i\), \(k - 1\)]]\)(...(\!\(\*SubscriptBox[\(s\), SubscriptBox[\(i\), \(1\)]]\)(k)).";
ReducedWord::usage="ReducedWord[w] returns a a reduced word {\!\(\*SubscriptBox[\(i\), \(1\)]\),...,\!\(\*SubscriptBox[\(i\), \(k\)]\)} such that for \!\(\*SubscriptBox[\(s\), \(j\)]\) the transposition (j j+1), w(k)=\!\(\*SubscriptBox[\(s\), SubscriptBox[\(i\), \(k\)]]\)(\!\(\*SubscriptBox[\(s\), SubscriptBox[\(i\), \(k - 1\)]]\)(...(\!\(\*SubscriptBox[\(s\), SubscriptBox[\(i\), \(1\)]]\)(k)).";
DividedDifference::usage ="DividedDifference[i,f] returns the ith divided difference of the polynomial f in variables {x[j]}.";
IsobaricDividedDifference::usage="IsobaricDividedDifference[i,f] returns the ith isobaric divided difference of the polynomial f in variables {x[j]}.";
LehmerCode::usage = "LehmerCode[w] returns the Lehmer code of the permutation w.";
LehmerCodeToPermutation::usage = "LehmerCodeToPermutation[c] returns the permutation w with Lehmer code c.";
OneFixedDominantQ::usage = "OneFixedDominantQ[w] returns True if w is of the form w=1w' where w' is dominant on {2,3,...,n} and False otherwise." ;
CoreRegion::usage = "CoreRegion[w] returns the core region of the permutation w, the possible squares where crosses can occur in pipe dreams.";
RotheDiagram::usage = "RotheDiagram[w] returns the Rothe diagram of the permutation w, the set of all (i,j) in [n\!\(\*SuperscriptBox[\(]\), \(2\)]\) such that w(i)>j and \!\(\*SuperscriptBox[\(w\), \(-1\)]\)(j)>i.";
BottomReducedPipeDream::usage = "BottomReducedPipeDream[w] returns the bottom reduced pipe dream of the permutation w.";
TopReducedPipeDream::usage = "TopReducedPipeDream[w] returns the top reduced pipe dream of the permutation w.";
AvoidsPattern::usage="AvoidsPattern[\[Sigma],\[Tau]] checks if the permutation sigma avoids the permutation pattern tau.";
GeneralizedPermutahedronZVector::usage="GeneralizedPermutahedronZVector[pointlist] returns the \!\(\*SubscriptBox[\(z\), \(I\)]\) values for the smallest generalized permutahedron containing pointlist.";
GeneralizedPermutahedronInequalities::usage="GeneralizedPermutahedronInequalities[pointlist] returns the defining inequalities for the smallest generalized permutahedron containing pointlist.";
GeneralizedPermutahedronYVector::usage="GeneralizedPermutahedronYVector[pointlist] returns the \!\(\*SubscriptBox[\(y\), \(I\)]\) values for the smallest generalized permutahedron containing pointlist.";
MConvexSetQ::usage="MConvexSetQ[J] checks whether or not J is an M-convex set, that is if J is exactly the set of integer points of a generalized permutahedron.";
Coefficients::usage="Coefficients[poly,vars] returns a list of the coefficients of the monomials occuring in the polynomial poly relative to the variables vars.\nCoefficients[poly] returns a list of the coefficients of monomials occuring in the polynomial poly relative to the variables detected in poly.";
Exponents::usage="Exponents[poly,vars] returns a list of the exponent vectors of the polynomial poly relative to the variables vars.\nExponents[poly] returns a list of the exponent vectors of the polynomial poly relative to the variables detected in poly.";
Monomials::usage="Monomials[poly,vars] returns a list of the terms of the polynomial poly relative to the variables vars.\nMonomials[poly] returns a list of the terms of the polynomial poly relative to the variables detected in poly.";
DominantQ::usage="DominantQ[w] checks whether the permutation w is dominant.";
GrassmannianQ::usage="GrassmannianQ[w] checks whether the permutation w is Grassmannian.";
VexillaryQ::usage="VexillaryQ[w] checks whether the permutation w is vexillary.";
DrawRotheDiagram::usage="DrawRotheDiagram[w] prints a graphical representation of the Rothe diagram of the permutation w.";
DemazureDifference::usage="DemazureDifference[i,f] returns the ith Demazure difference of the polynomial f in variables {x[j]}.";
DemazureLascouxDifference::usage="DemazureLascouxDifference[i,f] returns the ith Demazure-Lascoux difference of the polynomial f in variables {x[j]}.";
KeyPolynomial::usage="KeyPolynomial[alpha] returns the key polynomial of the composition alpha.";
HomogeneousComponents::usage="HomogeneousComponents[poly,vars] returns a list of the homogeneous components of the polynomial poly relative to the variables vars.\nHomogeneousComponents[poly] returns a list of the homogeneous components of the polynomial poly relative to the variables detected in poly.";
SymmetricPolynomialQ::usage="SymmetricPolynomialQ[poly,vars] returns true if the polynomial poly is symmetric with respect to the variables vars.";
RankMatrix::usage="RankMatrix[w] returns the rank matrix of the permutation w.";
BruhatOrderLessEqualQ::usage="BruhatOrderLessEqualQ[w,v] returns True if w\[LessEqual]v in the strong Bruhat order on \!\(\*SubscriptBox[\(S\), \(n\)]\).";
MultiplicityPolynomials::usage="MultiplicityPolynomials[poly,vars] returns the set of polynomials obtained by decrementing the coefficients of poly until they are all zero.";
XVariables::usage="XVariables[k] returns the list of variables {x[1],x[2],...,x[k]} for k a nonnegative integer.\nxvariables[w] returns the list of variables {x[1],x[2],...,x[n-1]} for w a permutation of [n].";
XBVariables::usage="XBVariables[k] returns the list of variables {x[1],x[2],...,x[k],b} for k a nonnegative integer.\nxbvariables[w] returns the list of variables {x[1],x[2],...,x[n-1],b} for w a permutation of [n].";
XYVariables::usage="XYVariables[k] returns the list of variables {x[1],x[2],...,x[k],y[1],y[2],...,y[k]} for k a nonnegative integer.\nxyvariables[w] returns the list of variables {x[1],x[2],...,x[n-1],y[1],y[2],...,y[k]} for w a permutation of [n].";
XYBVariables::usage="XYBVariables[k] returns the list of variables {x[1],x[2],...,x[k],y[1],y[2],...,y[k],b} for k a nonnegative integer.\nxybvariables[w] returns the list of variables {x[1],x[2],...,x[n-1],y[1],y[2],...,y[k],b} for w a permutation of [n].";
BruhatOrderLessEqualCoverQ::usage="bruhatLequalCover[w,v] returns True if w<v is a cover relation in the strong Bruhat order on \!\(\*SubscriptBox[\(S\), \(n\)]\).";
SkylineDiagram::usage="SkylineDiagram[alpha] returns the skyline diagram of the composition alpha.";
BalancedLabelings::usage="BalancedLabelings[D] returns the set of balanced labellings of the diagram D each viewed as an integer matrix.";
PartitionQ::usage="PartitionQ[\[Lambda]] checks that \[Lambda] is a decreasing sequence of nonnegative numbers.";
SchurPolynomial::usage="SchurPolynomial[\[Lambda], n] returns the Schur polynomial of the partition \[Lambda] in variables {x[1],...,x[n]}.\nSchurPolynomial[\[Lambda]] returns the Schur polynomial of the partition \[Lambda] in variables {x[1],...,x[m]}, where \[Lambda] has m parts.";
WiringDiagram::usage="WiringDiagram[p] returns the wiring diagram of the permutation word p=(\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(m\)]\)), where the \!\(\*SubscriptBox[\(p\), \(i\)]\) are positive integers.";
Schubitope::usage="For w in \!\(\*SubscriptBox[\(S\), \(n\)]\), Schubitope[w] returns the submodular function on subsets of [n] defining the Schubitope of w. ";
SchubitopeDimension::usage="SchubitopeDimension[w] returns the dimension of the Schubitope of w.";
ToZeroOneMatrix::usage="ToZeroOneMatrix[D] returns a matrix representation D, where D is a subset of [n] x [n].";
KohnertMoveResults::usage="KohnertMoveResults[D] returns the set of Kohnert diagrams generated by a Kohnert move on D, the matrix of a diagram.";
KohnertDiagrams::usage="KohnertDiagrams[D] returns the set of all Kohnert diagrams generated by any sequences of Kohnert moves on D, the matrix of a diagram.";
KohnertWeight::usage="KohnertWeight[D] returns the vector weight of D, the matrix of a diagram.";
KohnertPolynomial::usage="KohnertPolynomial[D] returns the polynomial generated by Kohnert diagrams of D, the matrix of a diagram.";
RandomDiagram::usage="RandomDiagram[m,n] returns a random subset of [m] x [n] as a zero-one matrix.";
DrawKohnertDiagram::usage="DrawKohnertDiagram[D] gives a visual representation of the Kohnert diagram D, repersentated by a matrix.";
DrawDiagram::usage="DrawDiagram[D, {m,n}] takes a subset D of [m] x [n] and represents it visually in a grid.";
ApplyKohnertMove::usage="ApplyKohnertMove[D,r] applies a Kohnert move to row r of the diagram D viewed as a matrix.";
ApplyKohnertSequence::usage="ApplyKohnertSequence[D,s] applies the sequences s of Kohnert moves to the diagram D viewed as a matrix.";
ValidKohnertMoveQ::usage="ValidKohnertMoveQ[D,r] checks that row r admits a valid Kohnert move in D.";
ValidKohnertSequenceQ::usage="ValidKohnertSequenceQ[D,s] checks that s is a valid sequence of Kohnert moves on D.";
ValidKohnertMoves::usage="ValidKohnertMoves[D] returns the list of all rows in D admitting a Kohnert move.";
ValidKohnertSequences::usage="ValidKohnertSequences[D] returns the list of all valid sequences of Kohnert moves on D.";
ValidRightmostLadderMoveQ::usage="ValidRightmostLadderMoveQ[P,r] checks that row r of the pipe dream P admits a ladder move.";
ApplyRightmostLadderMove::usage="ApplyRightmostLadderMove[P,r] applies the rightmost ladder move in row r to the pipe dream P.";
ApplyRightmostLadderMoveSequence::usage="ApplyRightmostLadderMoveSequence[P,s] applies the sequence s of rightmost ladder moves to the pipe dream P.";
ValidRightmostLadderMoves::usage="ValidRightmostLadderMoves[P] returns the list of rows of P admitting a rightmost ladder move.";
ValidRightmostLadderMoveSequenceQ::usage="ValidRightmostLadderMoveSequenceQ[P,s] checks that the sequence s of rightmost moves is a valid sequence on the pipe dream P.";
ValidRightmostLadderMoveSequences::usage="ValidRightmostLadderMoveSequences[P] returns the list of all valid rightmost ladder move sequences on the pipe dream P.";
PipeDreamWeight::usage="PipeDreamWeight[P] returns the weight of the pipe dream P, that is the vector whose \!\(\*SuperscriptBox[\(i\), \(th\)]\) entry is number of crosses in row i of P.";
Mitosis::usage="Mitosis[i,P] applies the \!\(\*SuperscriptBox[\(i\), \(th\)]\) mitosis algorithm to the pipe dream P. Setting Reduced->True or False specifices whether to return reduced pipe dreams or all pipe dreams.";
FultonEssentialSet::usage="FultonEssentialSet[w] returns the Fulton essential set of the permutation w, that is the boxes that are SE corners in the Rothe diagram of w.";
FundamentalSlidePolynomial::usage="FundamentalSlidePolynomial[\[Alpha]] returns the fundamental slide polynomial of alpha.";
Compositions::usage="Compositions[n,k] returns the integer compositions of n with at most k parts.\nCompositions[n,{k}] returns the integer compositions of n with exactly k parts.";
FundamentalSlideExpansion::usage="FundamentalSlideExpansion[w] returns the list of compositions whose fundamental slide polynomials occur in the expansion of the Schubert polynomial of w.";
MinimalPermutation::usage="MinimalPermutation[w] returns {w(1),...,w(n)} where w(i)=i for all i>n.";
KeyExpansion::usage="KeyExpansion[w] returns the set of compositions whose key polynomials sum to the Schubert polynomial of w.";
SchubertExpansion::usage="SchubertExpansion[P] returns the list of permutations and list of coefficients that give the polynomial P in variables x[i] as a linear combination of Schubert polynomials. Optional flag Memoization can be set to True or False to enable automatic storing of output.";
NorthWestDiagramMatrixQ::usage="NorthWestDiagramMatrixQ[M] checks if M is the 01 matrix corresponding to a northwest diagram.";
StronglySeparatedDiagramMatrixQ::usage="StronglySeparatedDiagramMatrixQ[M] checks if M is the 01 matrix corresponding to a strongly separated diagram.";
SchubertPolynomialQ::usage="SchubertPolynomialQ[poly] checks whether or not poly is a Schubert polynomial in variables x[i]. Optional flag Memoization can be set to True or False to enable automatic storing of output.";
PercentageAvoidingDiagramMatrixQ::usage="PercentageAvoidingDiagramMatrixQ[M] checks if M is the 01 matrix corresponding to a %-avoiding diagram.";
MConvexSupportQ::usage="MConvexSupportQ[f] tests whether the support of the polynomial f is M-convex. That is, whether the Newton polytope of f is a generalized permutahedron and is saturated by f.";
NormalizePolynomial::usage="NormalizePolynomial[f] returns the normalization of the polynomial f. Specifically, each term is divided by the product of factorials of its exponent vector.";
QuadraticFormToMatrix::usage="QuadraticFormToMatrix[f,n] takes a quadratic form f in variables x[i] and the number n of such variables as inputs and computes the n by n matrix corresponding to f.";
LorentzianPolynomialQ::usage="LorentzianPolynomialQ[h] checks whether the polynomial h is Lorentzian. That is, whether h has nonnegative coefficients, M-convex support, and satisfies the relevant derivative condition.";
CheckNonnegativity::usage="CheckNonnegativity is an option to LorentzianPolynomialQ, that can be set to True or False to enable or disable the checking of nonnegative coefficients in the input polynomial. The defualt value is True.";
CheckMConvexity::usage="CheckMConvexity is an option to LorentzianPolynomialQ, that can be set to True or False to enable or disable the checking of M-convexity of the support of the input polynomial. The defualt value is True.";
HeckeReduce::usage="HeckeReduce[r] reduces the integer word r using the nilHecke relations \!\(\*SuperscriptBox[SubscriptBox[\(s\), \(i\)], \(2\)]\)=\!\(\*SubscriptBox[\(s\), \(i\)]\), \!\(\*SubscriptBox[\(s\), \(i\)]\)\!\(\*SubscriptBox[\(s\), \(i + 1\)]\)\!\(\*SubscriptBox[\(s\), \(i\)]\)=\!\(\*SubscriptBox[\(s\), \(i + 1\)]\)\!\(\*SubscriptBox[\(s\), \(i\)]\)\!\(\*SubscriptBox[\(s\), \(i + 1\)]\), and \!\(\*SubscriptBox[\(s\), \(i\)]\)\!\(\*SubscriptBox[\(s\), \(j\)]\)=\!\(\*SubscriptBox[\(s\), \(j\)]\)\!\(\*SubscriptBox[\(s\), \(i\)]\) (if |i-j|>1).";
PolynomialDegree::usage="PolynomialDegree[f] returns the degree of the multivariable polynomial f.";
RotheBPD::usage="RotheBPD[w] returns Rothe bumpless pipe dream for a given permutation w.";
ASMQ::usage="ASMQ[P] tests whether a given square matrix P is an alternating sign matrix.";
DroopMove::usage="DroopMove[P,i0,j0,imax,jmax] returns a pipe dream which differs from P by a droop move performed in the rectangle defined by the corners (i0,j0) and (imax, jmax).";
DroopMoveQ::usage="DroopMoveQ[P,i0,j0,imax,jmax] tests whether it is possible to perform a valid droop move on the pipe dream P within the rectangle defined by the corners (i0,j0) and (imax, jmax).";
ReducedBumplessPipeDreams::usage="ReducedBumplessPipeDreams[w] generates all the reduced bumpless pipe dreams for the permutation w.";
BPDtoPermutation::usage="BPDtoPermutation[P] returns the permutation associated with the bumpless pipe dream P.";
DrawBumplessPipeDream::usage="DrawBumplessPipeDream[P] generates a visual representation of the bumpless pipe dream P (given as an alternating sign matrix).";
PipeColor::usage="PipeColor is an option to DrawBumplessPipeDream and DrawMarkedBumplessPipeDream that can be set to any color.  The default value is green.";
kDroop1::usage="kDroop1[P,i0,j0,imax,jmax] returns a bumpless pipe dream which differs from P by a k-theoretic droop move performed in the rectangle defined by the corners (i0,j0) and (imax, jmax).";
KDroop1Q::usage="KDroop1Q[P,i0,j0,imax,jmax] tests whether it is possible to perform a k-theoretic droop move on the bumpless pipe dream P in the rectangle defined by the corners (i0,j0) and (imax, jmax).";
kDroop2::usage="kDroop2[P,i0,j0,imax,jmax] returns a bumpless pipe dream which differs from P by a k-theoretic droop move performed in the rectangle defined by the corners (i0,j0) and (imax, jmax).";
KDroop2Q::usage="KDroop2Q[P,i0,j0,imax,jmax] tests whether it is possible to perform a k-theoretic droop move on the bumpless pipe dream P in the rectangle defined by the corners (i0,j0) and (imax, jmax).";
BumplessPipeDreams::usage="BumplessPipeDreams[w] returns all bumpless pipe dreams (reduced and nonreduced) for the permutation w.";
PivotQ::usage="PivotQ[P,a,b,i,j] determines whether the entry (i,j) is a pivot of (a,b) in the alternating sign matrix P.";
RemovableNeg1Q::usage="RemovableNeg1Q[P,a,b] determines whether the entry (a,b) in the alternating sign matrix P is a removable -1.";
BPDtoDemazureProduct::usage="BPDtoDemazureProduct[P] outputs the permutation correponding to the Demazure product of the column reading word of the bumpless pipe dream P.";
MarkedBumplessPipeDreams::usage="MarkedBumplessPipeDreams[w] returns all marked bumpless pipe dreams for the permutation w.";
DrawMarkedBumplessPipeDream::usage="DrawMarkedBumplessPipeDream[{P,s}] generates a visual representation of the input marked bumpless pipe dream.";

Unprotect[
Memoization,
PermutationListForm,
Permutations,
PipeDreamQ,
PipeDreamToPermutation,
PipeDreams,
DrawPipeDream,
Inversions,
NumberInversions,
ReducedPipeDreamQ,
ReducedPipeDreams,
GrothendieckPolynomial,
BetaGrothendieckPolynomial,
DoubleGrothendieckPolynomial,
BetaDoubleGrothendieckPolynomial,
SchubertPolynomial,
Descents,
Ascents,
Transposition,
TransitionRule,
DoubleSchubertPolynomial,
Letter,
ReducedWords,
ReducedWord,
DividedDifference,
IsobaricDividedDifference,
LehmerCode,
LehmerCodeToPermutation,
OneFixedDominantQ,
CoreRegion,
RotheDiagram,
BottomReducedPipeDream,
TopReducedPipeDream,
AvoidsPattern,
GeneralizedPermutahedronZVector,
GeneralizedPermutahedronInequalities,
GeneralizedPermutahedronYVector,
MConvexSetQ,
Coefficients,
Exponents,
DominantQ,
GrassmannianQ,
VexillaryQ,
DrawRotheDiagram,
DemazureDifference,
DemazureLascouxDifference,
KeyPolynomial,
HomogeneousComponents,
SymmetricPolynomialQ,
RankMatrix,
BruhatOrderLessEqualQ,
MultiplicityPolynomials,
XVariables,
XBVariables,
XYVariables,
XYBVariables,
BruhatOrderLessEqualCoverQ,
SkylineDiagram,
BalancedLabelings,
PartitionQ,
SchurPolynomial,
WiringDiagram,
Schubitope,
SchubitopeDimension,
ToZeroOneMatrix,
KohnertMoveResults,
KohnertDiagrams,
KohnertWeight,
KohnertPolynomial,
RandomDiagram,
DrawKohnertDiagram,
GhostKohnertMoveResults,
GhostKohnertDiagrams,
GhostKohnertWeight,
GhostKohnertPolynomial,
DrawDiagram,
ApplyKohnertMove,
ApplyKohnertSequence,
ValidKohnertMoveQ,
ValidKohnertSequenceQ,
ValidKohnertMoves,
ValidKohnertSequences,
ValidRightmostLadderMoveQ,
ApplyRightmostLadderMove,
ApplyRightmostLadderMoveSequence,
ValidRightmostLadderMoves,
ValidRightmostLadderMoveSequenceQ,
ValidRightmostLadderMoveSequences,
PipeDreamWeight,
Mitosis,
FultonEssentialSet,
FundamentalSlidePolynomial,
Compositions,
FundamentalSlideExpansion,
MinimalPermutation,
KeyExpansion,
SchubertExpansion,
NorthWestDiagramMatrixQ,
StronglySeparatedDiagramMatrixQ,
SchubertPolynomialQ,
PercentageAvoidingDiagramMatrixQ,
SkewSchurPolynomial,
MConvexSupportQ,
NormalizePolynomial,
QuadraticFormToMatrix,
LorentzianPolynomialQ,
CheckNonnegativity,
CheckMConvexity,
HeckeReduce,
PolynomialDegree,
RotheBPD,
ASMQ,
DroopMove,
DroopMoveQ,
ReducedBumplessPipeDreams,
BPDtoPermutation,
DrawBumplessPipeDream,
PipeColor,
kDroop1,
KDroop1Q,
kDroop2,
KDroop2Q,
BumplessPipeDreams,
PivotQ,
RemovableNeg1Q,
BPDtoDemazureProduct,
MarkedBumplessPipeDreams,
DrawMarkedBumplessPipeDream
];

Unprotect[Global`x,Global`y,Global`b,Global`t];
Global`x=Global`x;
Global`y=Global`y;
Global`b=Global`b;
Global`t=Global`t;
Protect[Global`x,Global`y,Global`b,Global`t];

Begin["`Private`"];
x=Global`x;
y=Global`y;
b=Global`b;
t=Global`t;

BoxListQ[boxlist_]:=Module[{},
If[!ListQ[boxlist],Return[False]];
If[!AllTrue[boxlist,ListQ[#]&&Length[#]==2&&AllTrue[#,IntegerQ[#]&&#>0&]&],
Return[False];
];
Return[True];
];

PermutationListForm[w_?PermutationCyclesQ]:=Module[{},
If[w==Cycles[{}],Return[{1}]];
Return[PermutationList[w]];
];

Permutations[n_Integer?Positive]:=Return[System`Permutations[Range[n]]];

PipeDreamQ[P_?MatrixQ]:=Module[{zeroone,sezeros},
zeroone=MatrixQ[P,#==1||#==0&];
sezeros=SubsetQ[Position[P,0],Select[Tuples[Range[Length[P]],2],#[[1]]+#[[2]]>=Length[P]+1&]];
Return[sezeros&&zeroone];
];

PipeDreamToPermutation[P_?PipeDreamQ]:=Module[{info,crossings,i,j,n},
n=Length[P];
info=Table[{0,0},{n},{n}];
crossings={};
For[i=n,i>=1,i--,
For[j=1,j<=n-i+1,j++,
If[{i,j}=={n,1},info[[n,1]]={n,0}; ];
If[i==n&&j>1,info[[i,j]]=0];
If[j==1&&i<n&& P[[i,j]]==0,info[[i,j]]={i,info[[i+1,j]][[1]]};];
If[j==1&&i<n&& P[[i,j]]==1,info[[i,j]]={info[[i+1,j]][[1]],i}; AppendTo[crossings,Sort[info[[i,j]]]]];
If[i<n&&j>1&&P[[i,j]]==0,info[[i,j]]={info[[i,j-1]][[2]],info[[i+1,j]][[1]]}];
If[i<n&&j>1&&P[[i,j]]==1&&MemberQ[crossings,Sort[{info[[i,j-1]][[2]],info[[i+1,j]][[1]]}]],info[[i,j]]={info[[i,j-1]][[2]],info[[i+1,j]][[1]]}];
If[i<n&&j>1&&P[[i,j]]==1&&!MemberQ[crossings,Sort[{info[[i,j-1]][[2]],info[[i+1,j]][[1]]}]],info[[i,j]]={info[[i+1,j]][[1]],info[[i,j-1]][[2]]}; AppendTo[crossings,Sort[info[[i,j]]]]]
];
];
Return[InversePermutation[info[[1,All,1]]]];
];

(*
PipeDreams[w_?PermutationListQ]:=Module[{changeoneelbow,redpipes,newpipedreams,allpipes},
changeoneelbow[P_]:=Module[{elbows,Q},
elbows=Select[Position[P,0],Total[#]<=Length[P]&];
Return[Table[Q=P;Q[[elb[[1]],elb[[2]]]]=1;Q,{elb,elbows}]];
];
redpipes=ReducedPipeDreams[w];
newpipedreams=Flatten[changeoneelbow/@redpipes,1];
newpipedreams=Select[newpipedreams,PipeDreamToPermutation[#]==w&];
allpipes=DeleteDuplicates[Join[redpipes,newpipedreams]];
While[Length[newpipedreams]>0,
newpipedreams=Flatten[changeoneelbow/@newpipedreams,1];
newpipedreams=Select[newpipedreams,PipeDreamToPermutation[#]==w&];
allpipes=DeleteDuplicates[Join[allpipes,newpipedreams]];
];
Return[Sort[allpipes]];
];
*)

PipeDreams[w_?PermutationListQ]:=Module[{w0,P0,word,pd,v},
w0=Reverse[Range[Max[w]]];
P0=BottomReducedPipeDream[w0];
word=ReducedWord[PermutationProduct[w0,InversePermutation[w]]];
pd={P0};
v=w0;
Do[
v=MinimalPermutation[PermutationProduct[Letter[i],v]];
pd=Mitosis[i,#,Reduced->False]&/@pd;
pd=DeleteCases[pd,{}];
pd=Flatten[Select[pd,MinimalPermutation[PipeDreamToPermutation[#[[1]]]]==v&],1];
,{i,word}];
Return[pd];
]

Options[DrawPipeDream]={CrossColor->Purple,ElbowColor->Purple,Labels->True};
DrawPipeDream[P_?PipeDreamQ,OptionsPattern[]]:=Module[{vertgridlines,horizgridlines,Pcrosspositions,changeofcoordinates,Crossing,crosses,Pelbowpositions,normalelbowpositions,boundaryelbowposition,NormalElbow,normalelbows,BoundaryElbow,boundaryelbows,elbows,upperlabels,perm,leftlabels,picture},
vertgridlines=Table[{Thick,Line[{{0,j},{Length[P],j}}]},{j,0,Length[P]}];
horizgridlines=Table[{Thick,Line[{{j,0},{j,Length[P]}}]},{j,0,Length[P]}];
Pcrosspositions=Position[P,1];
changeofcoordinates[{i_,j_}]:={j-1/2,Length[P]+1/2-i};
If[Pcrosspositions!={},
Crossing[{i_,j_}]:={{Thick,OptionValue[CrossColor],Line[{{i-1/2,j},{i+1/2,j}}]},{Thick,OptionValue[CrossColor],Line[{{i,j-1/2},{i,j+1/2}}]}};
crosses=Crossing/@(changeofcoordinates/@Pcrosspositions);
,crosses={};
];
Pelbowpositions=Position[P,0];
normalelbowpositions=Select[Pelbowpositions,#[[1]]+#[[2]]<Length[P]+1&];
boundaryelbowposition=Select[Pelbowpositions,#[[1]]+#[[2]]==Length[P]+1&];
If[normalelbowpositions!={},
NormalElbow[{i_,j_}]:={{Thick,OptionValue[ElbowColor],Circle[{i+1/2,j-1/2},1/2,{Pi/2,Pi}]},{Thick,OptionValue[ElbowColor],Circle[{i-1/2,j+1/2},1/2,{0,-Pi/2}]}};
normalelbows=NormalElbow/@(changeofcoordinates/@normalelbowpositions);
,normalelbows={};
];
BoundaryElbow[{i_,j_}]:={{Thick,OptionValue[ElbowColor],Circle[{i-1/2,j+1/2},1/2,{0,-Pi/2}]}};
boundaryelbows=BoundaryElbow/@(changeofcoordinates/@boundaryelbowposition);
elbows=Join[normalelbows,boundaryelbows];
upperlabels=Table[{Text[Style[i,Large],{i-1/2,Length[P]+1/5}]},{i,1,Length[P]}];
perm=PipeDreamToPermutation[P];
leftlabels=Table[{Text[Style[perm[[i]],Large],{-1/5,Length[P]+1/2-i}]},{i,1,Length[P]}];
If[OptionValue[Labels],
picture=Graphics[Join[vertgridlines,horizgridlines,Flatten[crosses,1],Flatten[elbows,1],upperlabels,leftlabels]],
picture=Graphics[Join[vertgridlines,horizgridlines,Flatten[crosses,1],Flatten[elbows,1]]]
];
Return[picture];
];

Inversions[w_?PermutationListQ]:=Return[Select[Tuples[Range[Length[w]],{2}],#[[1]]<#[[2]]&&w[[#[[1]]]]>w[[#[[2]]]]&]];

NumberInversions[w_?PermutationListQ]:=Module[{pairs},
pairs=Select[Tuples[Range[Length[w]],{2}],#[[1]]<#[[2]]&];
Return[Length[Select[pairs,w[[#[[1]]]]>w[[#[[2]]]]&]]];
];

ReducedPipeDreamQ[P_?PipeDreamQ]:=Return[Length[Position[P,1]]==NumberInversions[PipeDreamToPermutation[P]]]

ReducedPipeDreams[w_?PermutationListQ]:=Module[{w0,P0,word,rpd},
w0=Reverse[Range[Max[w]]];
P0=BottomReducedPipeDream[w0];
word=ReducedWord[PermutationProduct[w0,InversePermutation[w]]];
rpd={P0};
Do[
rpd=Flatten[Mitosis[i,#,Reduced->True]&/@rpd,1];
,{i,word}];
Return[rpd];
];

Letter[i_Integer?Positive]:=Return[Join[Table[j,{j,1,i-1}],{i+1,i}]];

ReducedWord[w_?PermutationListQ]:=Module[{descents,u,word},
If[Sort[w]==w,Return[{}]];
descents=Select[Range[Length[w]-1],w[[#]]>w[[#+1]]&];
u=PermutationProduct[Letter[descents[[1]]],w];
Return[Prepend[ReducedWord[u],descents[[1]]]];
];

ReducedWords[w_?PermutationListQ]:=Module[{descents,newperms,newwords},
If[Sort[w]==w,Return[{{}}]];
descents=Select[Range[Length[w]-1],w[[#]]>w[[#+1]]&];
newperms=Table[PermutationProduct[Letter[i],w],{i,descents}];
newwords=ReducedWords/@newperms;
Return[Flatten[Table[(Prepend[#,descents[[i]]])&/@newwords[[i]],{i,1,Length[descents]}],1]];
];

DividedDifference[i_Integer?Positive,f_]:=Module[{diff,vars},
vars=Variables[f];
diff=PolynomialReduce[(f-(f/.{x[i]->x[i+1],x[i+1]->x[i]})),{(x[i]-x[i+1])},vars];
Return[Expand[diff[[1,1]]]];
];

IsobaricDividedDifference[i_Integer?Positive,f_]:=Return[Expand[DividedDifference[i,(1-x[i+1])f]]];

Descents[w_?PermutationListQ]:=Return[Select[Range[Length[w]-1],w[[#]]>w[[#+1]]&]];

Ascents[w_?PermutationListQ]:=Return[Select[Range[Length[w]-1],w[[#]]<w[[#+1]]&]];

Transposition[i_Integer?Positive,j_Integer?Positive]:=(Range[Max[{i,j}]]/.{i->j,j->i});

TransitionRule[w_?PermutationListQ]:=Module[{r,s,v,transitions,transperms},
If[w==Range[Length[w]],Return[{}]];
r=Max[Descents[w]];
s=Max[Select[Range[r+1,Length[w]],w[[#]]<w[[r]]&]];
v=PermutationProduct[Transposition[r,s],w];
transitions=Select[Range[r-1],NumberInversions[PermutationProduct[Transposition[#,r],v]]==NumberInversions[w]&];
transperms=Table[PermutationProduct[Transposition[q,r],v],{q,transitions}];
Return[{{r,v},transperms}];
];

BetaGrothendieckPolynomial[w_?PermutationListQ]:=Module[{isobaricDividedDifferenceBeta,w0,u,word,g,i},
isobaricDividedDifferenceBeta[i_Integer?Positive,f_]:=Return[Expand[DividedDifference[i,(1+b x[i+1])f]]];
w0=Reverse[Range[Length[w]]];
u=PermutationProduct[w0,InversePermutation[w]];word=ReducedWord[u];
g=Apply[Times,Array[x[#]^(Length[w]-#)&,Length[w]-1]];
For[i=1,i<=Length[word],i++,g=isobaricDividedDifferenceBeta[word[[i]],g];
];
Return[Expand[g]];
];

BetaDoubleGrothendieckPolynomial[w_?PermutationListQ]:=Module[{isobaricDividedDifferenceBeta,w0,u,word,g,i},
isobaricDividedDifferenceBeta[i_Integer?Positive,f_]:=Return[Expand[DividedDifference[i,(1+b x[i+1])f]]];
w0=Reverse[Range[Length[w]]];
u=PermutationProduct[w0,InversePermutation[w]];
word=ReducedWord[u];
g=Apply[Times,(x[#[[1]]]+y[#[[2]]]-x[#[[1]]]y[#[[2]]])&/@Select[Tuples[Range[Length[w]-1],{2}],#[[1]]+#[[2]]<=Length[w]&]];
For[i=1,i<=Length[word],i++,
g=isobaricDividedDifferenceBeta[word[[i]],g];
];
Return[Expand[g]];
];

LehmerCode[w_?PermutationListQ]:=Return[Table[Count[w[[#]]<w[[i]]&/@Range[i+1,Length[w]],True],{i,1,Length[w]}]];

LehmerCodeToPermutation[list_?CompositionQ]:=Module[{k,begin,l,i,j},
k=list;
Label[begin];
l=k;
For[i=Length[l]-1,i>=1,i--,
For[j=i+1,j<=Length[l],j++,
If[l[[j]]>=l[[i]], l[[j]]=l[[j]]+1;];
];
];
If[PermutationListQ[l+1],Return[l+1]];
AppendTo[k,0];
Goto[begin];
];

OneFixedDominantQ[w_?PermutationListQ]:=Module[{code,decreasing},
If[w[[1]]!=1, Return[False]];
code=LehmerCode[w];
decreasing=Table[code[[i]]>=code[[i+1]],{i,2,Length[code]-1}];
If[DeleteDuplicates[decreasing]!={True},Return[False]];
Return[True];
];

CoreRegion[w_?PermutationListQ]:=Module[{pipes},
pipes=PipeDreams[w];
Return[DeleteDuplicates[Flatten[Table[Position[P,1],{P,pipes}],1]]];
];

RotheDiagram[w_?PermutationListQ]:=Module[{allSquares},
allSquares=Tuples[Range[Length[w]],{2}];
Return[Select[allSquares,w[[#[[1]]]]>#[[2]]&&InversePermutation[w][[#[[2]]]]>#[[1]]&]];
];

BottomReducedPipeDream[w_?PermutationListQ]:=Module[{rothe,perRow,P,i,j},
rothe=RotheDiagram[w];
perRow=Table[Length[Select[rothe,#[[1]]==j&]],{j,1,Length[w]}];
P=Table[0,{Length[w]},{Length[w]}];
For[j=1,j<=Length[perRow],j++,
For[i=1,i<=perRow[[j]],i++,
P[[j,i]]=1;
];
];
Return[P];
];

TopReducedPipeDream[w_?PermutationListQ]:=Module[{rothe,perColumn,P,i,j},
rothe=RotheDiagram[w];
perColumn=Table[Length[Select[rothe,#[[2]]==j&]],{j,1,Length[w]}];
P=Table[0,{Length[w]},{Length[w]}];
For[j=1,j<=Length[perColumn],j++,
For[i=1,i<=perColumn[[j]],i++,
P[[i,j]]=1;
];
];
Return[P];
];

AvoidsPattern[sigma_?PermutationListQ,tau_?PermutationListQ]:=Module[{n,k,possibles,indices,pairs,list1,list2,ind,results},
n=Length[sigma];
k=Length[tau];
possibles=Subsets[Range[n],{k}];
indices=Tuples[Range[k],2];
pairs=Tuples[Range[k],{2}];
results=Table[
list1=Select[pairs,sigma[[current[[#[[1]]]]]]<sigma[[current[[#[[2]]]]]]&];
list2=Select[pairs,tau[[#[[1]]]]<tau[[#[[2]]]]&];
Sort[list1]==Sort[list2],{current,possibles}];
If[MemberQ[results,True],Return[False],Return[True]];
];

(*
GeneralizedPermutahedronInequalities[list_?MatrixQ]:=Module[{subsets,lowerbound,zvector,zz,inequalities},
subsets=Subsets[Range[Length[list[[1]]]]];
lowerbound[I_]:=Min[Map[Total,list[[All,I]]]];
zvector=lowerbound/@subsets;
zz[set_]:=zvector[[Position[subsets,set][[1,1]]]];
inequalities=Table[Sum[t[i],{i,I}]>=zz[I],{I,Delete[Delete[subsets,1],-1]}];
AppendTo[inequalities,Sum[t[i],{i,1,Length[list[[1]]]}]==zz[subsets[[-1]]]];
Return[inequalities];
];
*)

GeneralizedPermutahedronYVector[list_?MatrixQ]:=Module[{subsets,upperbound,zvector,zz,inequalities,answers,yvector},
subsets=Subsets[Range[Length[list[[1]]]]];
upperbound[I_]:=Max[Map[Total,list[[All,I]]]];
zvector=upperbound/@subsets;
zz[set_]:=zvector[[Position[subsets,set][[1,1]]]];
inequalities=Table[Sum[t[i],{i,I}]<=zz[I],{I,Delete[Delete[subsets,1],-1]}];
AppendTo[inequalities,Sum[t[i],{i,1,Length[list[[1]]]}]==zz[subsets[[-1]]]];
yvector=Table[Sum[(-1)^(Length[Complement[I,J]])zz[J],{J,Select[subsets,SubsetQ[I,#]&]}],{I,subsets}];
Return[yvector];
];

GeneralizedPermutahedronZVector[list_?MatrixQ]:=Module[{subsets,upperbound,zvector},
subsets=Subsets[Range[Length[list[[1]]]]];
upperbound[I_]:=Max[Map[Total,list[[All,I]]]];
zvector=upperbound/@subsets;
Return[zvector];
];

SubmodularVectorQ[zvect_?VectorQ]:=Module[{len,I,J,posI,posJ,posIuJ,posIJ,results,subsets,pairs,ddresults},
len=Log[2,Length[zvect]];
subsets=Subsets[Range[len]];
pairs=Tuples[subsets,{2}];
results=Reap[Do[{I,J}=set;
posI=Position[subsets,I][[1,1]];
posJ=Position[subsets,J][[1,1]];
posIuJ=Position[subsets,Union[I,J]][[1,1]];
posIJ=Position[subsets,Intersection[I,J]][[1,1]];
Sow[zvect[[posI]]+zvect[[posJ]]>=zvect[[posIJ]]+zvect[[posIuJ]]];,{set,pairs}]][[2,1]];
ddresults=DeleteDuplicates[results];
Return[Length[ddresults]==1&&ddresults[[1]]==True];
];

GeneralizedPermutahedronInequalities[list_?MatrixQ]:=Module[{subsets,upperbound,zvector,zz,inequalities},
subsets=Subsets[Range[Length[list[[1]]]]];
upperbound[I_]:=Max[Map[Total,list[[All,I]]]];
zvector=upperbound/@subsets;
zz[set_]:=zvector[[Position[subsets,set][[1,1]]]];
inequalities=Table[Sum[t[i],{i,I}]<=zz[I],{I,Delete[Delete[subsets,1],-1]}];
AppendTo[inequalities,Sum[t[i],{i,1,Length[list[[1]]]}]==zz[subsets[[-1]]]];
Return[inequalities];
];

MConvexSetQ[list_?(MatrixQ[#,IntegerQ]&)]:=Module[{subsets,upperbound,zvector,zz,inequalities,answers,saturated,submodular},
If[list=={},Return[True]];
If[!MatrixQ[list,NumberQ],Return[False]];
If[Length[DeleteDuplicates[Total/@list]]!=1,Return[False]];
subsets=Subsets[Range[Length[list[[1]]]]];
upperbound[I_]:=Max[Map[Total,list[[All,I]]]];
zvector=upperbound/@subsets;
submodular=SubmodularVectorQ[zvector];
If[!submodular,Return[False]];
zz[set_]:=zvector[[Position[subsets,set][[1,1]]]];
inequalities=Table[Sum[t[i],{i,I}]<=zz[I],{I,Delete[Delete[subsets,1],-1]}];
AppendTo[inequalities,Sum[t[i],{i,1,Length[list[[1]]]}]==zz[subsets[[-1]]]];
answers=Table[t[i],{i,1,Length[list[[1]]]}]/.Solve[inequalities,Table[t[i],{i,1,Length[list[[1]]]}],Integers];
saturated=(Sort[answers]==Sort[DeleteDuplicates[list]]);
Return[saturated];
];

Monomials[poly_,vars_?VectorQ]:=MonomialList[poly,vars];
Monomials[poly_]:=Monomials[poly,Variables[poly]];

Coefficients[poly_,vars_?VectorQ]:=Module[{replace2},
replace2[a_->b_]:=b;
Return[replace2/@CoefficientRules[poly,vars]];
];
Coefficients[poly_]:=Coefficients[poly,Variables[poly]];

Exponents[poly_,vars_?VectorQ]:=Module[{replace1},
replace1[a_->b_]:=a;
Return[replace1/@CoefficientRules[poly,vars]];
];
Exponents[poly_]:=Exponents[poly,Variables[poly]];

VexillaryQ[w_?PermutationListQ]:=Module[{rothe,columns},
rothe=RotheDiagram[w];
columns=DeleteCases[Table[Select[rothe,#[[2]]==j&],{j,1,Length[w]}],{}];
columns=#[[All,1]]&/@columns;
columns=Sort[columns,SubsetQ];
Return[AllTrue[Range[Length[columns]-1],SubsetQ[columns[[#]],columns[[#+1]]]&]];
];

GrassmannianQ[w_?PermutationListQ]:=Module[{descents},
descents=Select[Range[Length[w]-1],w[[#]]>w[[#+1]]&];
Return[Length[descents]==1];
];

DominantQ[w_?PermutationListQ]:=Module[{lc},
lc=LehmerCode[w];
Return[AllTrue[Range[Length[lc]-1],lc[[#]]>=lc[[#+1]]&]];
];

Options[DrawRotheDiagram]={Label->""};
DrawRotheDiagram[w_?PermutationListQ,OptionsPattern[]]:=
Module[{dots,vertgridlines,horizgridlines,horizlines,vertlines,squares,rothediagram,label},
dots=Table[{w[[i]]-1/2,Length[w]-i+1/2},{i,1,Length[w]}];
vertgridlines=Table[{Thick,Line[{{0,j},{Length[w],j}}]},{j,0,Length[w]}];
horizgridlines=Table[{Thick,Line[{{j,0},{j,Length[w]}}]},{j,0,Length[w]}];
horizlines={Thick,Red,#}&/@Table[Line[{dots[[i]],{dots[[i,1]],0}}],{i,1,Length[w]}];
vertlines={Thick,Red,#}&/@Table[Line[{dots[[i]],{Length[w],dots[[i,2]]}}],{i,1,Length[w]}];
squares={Purple,Rectangle[#[[1]],#[[2]]]}&/@({{#[[1]]-1,#[[2]]-1},{#[[1]],#[[2]]}}&/@({#[[2]],Length[w]-#[[1]]+1}&/@RotheDiagram[w]));
If[OptionValue[Label]=="",label={},label={Text[Style[OptionValue[Label],Large],{Length[w]/2,Length[w]+1}]}];
rothediagram=Graphics[Join[squares,horizgridlines,vertgridlines,({Red,PointSize[Large],#}&/@Point/@dots),{Point[{1,0}],horizlines,vertlines,label}]];
Return[rothediagram];
];

DemazureDifference[i_Integer?Positive,f_]:=Return[DividedDifference[i,x[i]f]];

DemazureLascouxDifference[i_Integer?Positive,f_]:=Return[DividedDifference[i,Expand[x[i](1-x[i+1])f]]]

CompositionQ[alpha_]:=Module[{},
Return[VectorQ[alpha,IntegerQ[#]&&#>=0&]];
];

KeyPolynomial[alpha_?CompositionQ]:=Module[{beta,descent,betahat},
beta=alpha;
While[Length[beta]<Max[beta],AppendTo[beta,0];];
If[Sort[beta,#1>#2&]==beta, Return[Apply[Times,Table[x[i],{i,1,Length[beta]}]^beta]]];
descent=Select[Range[Length[beta]-1],beta[[#+1]]>beta[[#]]&][[1]];
betahat=beta;
betahat[[descent]]=beta[[descent+1]];
betahat[[descent+1]]=beta[[descent]];
Return[DemazureDifference[descent,KeyPolynomial[betahat]]];
];

HomogeneousComponents[poly_,variables_?VectorQ]:=Module[{split,mults},
split=GatherBy[MonomialList[poly,variables],Total[Exponents[#,variables][[1]]]&];
mults=Flatten[(Total/@Exponents[#,variables])&/@split[[All,1]]];
split=Table[split[[Position[mults,i][[1,1]]]],{i,Sort[mults]}];
Return[Total/@split];
];
HomogeneousComponents[poly_]:=HomogeneousComponents[poly,Variables[poly]];

SymmetricPolynomialQ[polynomial_,variables_?VectorQ]:=Module[{swaps,tests},
swaps=Table[{variables[[i]]->variables[[i+1]],variables[[i+1]]->variables[[i]]},{i,1,Length[variables]-1}];
Return[AllTrue[Table[Exponents[(polynomial/.swap)-polynomial,variables]=={},{swap,swaps}],#&]];
];

RankMatrix[w_?PermutationListQ]:=Module[{r,indices,i,j},
r=ConstantArray[0,{Length[w],Length[w]}];
indices=Tuples[Range[Length[w]],{2}];
Do[
{i,j}=index;
r[[i,j]]=Length[Select[Range[i],w[[#]]<=j&]];
,{index,indices}];
Return[r];
];

BruhatOrderLessEqualQ[w_?PermutationListQ,v_?PermutationListQ]:=Module[{m,n,v1,w1},
m=Length[w];
n=Length[v];
If[m<n,
v1=v;
w1=Join[w,Table[i,{i,m+1,n}]];
];
If[m>n,
w1=w;
v1=Join[v,Table[i,{i,n+1,m}]];
];
If[m==n,
v1=v;
w1=w;
];
Return[AllTrue[Thread[Flatten[RankMatrix[v1]]<=Flatten[RankMatrix[w1]]],#&]];
];

MultiplicityPolynomials[poly_,vars_?VectorQ]:=Module[{exps,coeffs,numbershifts,shifts,shiftedcoeffs,i},
exps=Exponents[poly,vars];
coeffs=Coefficients[poly,vars];
numbershifts=Max[Abs/@coeffs];
shifts={Total[Table[coeffs[[i]]Apply[Times,vars^(exps[[i]])],{i,1,Length[exps]}]]};
shiftedcoeffs=coeffs;
For[i=1,i<numbershifts,i++,
shiftedcoeffs=shiftedcoeffs-Table[Piecewise[{{1,i>0},{-1,i<0},{0,i==0}}],{i,shiftedcoeffs}];
AppendTo[shifts,Total[Table[shiftedcoeffs[[i]]Apply[Times,vars^(exps[[i]])],{i,1,Length[exps]}]]];
];
Return[shifts];
];
MultiplicityPolynomials[poly_]:=MultiplicityPolynomials[poly,Variables[poly]];

XVariables[w_?PermutationListQ]:=Return[Array[x,Length[w]-1]];

XVariables[n_Integer?Positive]:=Return[Array[x,n]];

XBVariables[w_?PermutationListQ]:=Return[Append[Array[x,Length[w]-1],b]];

XBVariables[n_Integer?Positive]:=Return[Append[Array[x,n],b]];

XYVariables[w_?PermutationListQ]:=Return[Join[Array[x,Length[w]-1],Array[y,Length[w]-1]]];

XYVariables[n_Integer?Positive]:=Return[Join[Array[x,n],Array[y,n]]];

XYBVariables[w_?PermutationListQ]:=Return[Append[Join[Array[x,Length[w]-1],Array[y,Length[w]-1]],b]];

XYBVariables[n_Integer?Positive]:=Return[Append[Join[Array[x,n],Array[y,n]],b]];

BruhatOrderLessEqualCoverQ[w1_?PermutationListQ,v1_?PermutationListQ]:=Module[{v,w,n,diff,e,f},
n=Max[Join[v1,w1]];
w=Join[w1,Range[Max[w1]+1,n]];
v=Join[v1,Range[Max[v1]+1,n]];
If[NumberInversions[w]>=NumberInversions[v],Return[False]];
diff=Select[Range[n],w[[#]]!=v[[#]]&];
e=Min[diff];
f=Max[diff];
Return[w[[e]]<w[[f]]&&Length[Select[Range[e+1,f-1],w[[e]]<w[[#]]<w[[f]]&]]==0];
];

SkylineDiagram[alpha_?VectorQ]:=Flatten[DeleteCases[Table[Table[{j,i},{i,1,alpha[[j]]}],{j,1,Length[alpha]}],{}],1];

BalancedLabelings[diagram_?BoxListQ,n_Integer?Positive]:=Module[{hook,allhooks,balancedQ,fillingbalancedQ,variables,flagconditions,columns,strictcolumns,unbalancedfillings,balancedlabelings,A,matrices},
hook[square_]:=Select[diagram,#[[1]]==(square[[1]]&&#[[2]]>=square[[2]])||#[[2]]==square[[2]]&&#[[1]]>=square[[1]]&];
allhooks=hook/@diagram;
balancedQ[hook_,hookfilling_]:=Module[{base,right,down,goodperms,hookrightfilling,permincreasingright,hookdownfilling,permincreasingdown,reversedhook,sortedfilling},
base={Min[hook[[All,1]]],Min[hook[[All,2]]]};
right=Select[hook,#[[2]]>base[[2]]||#==base&];
down=Select[hook,#[[1]]>base[[1]]||#==base&];
reversedhook=Join[Reverse[Delete[right,1]],down];
sortedfilling=Sort[hookfilling];
Return[sortedfilling[[Position[reversedhook,base][[1,1]]]]==hookfilling[[Position[hook,base][[1,1]]]]];
];
fillingbalancedQ[filling_]:=Module[{inducedfilling},
inducedfilling[hook_]:=filling[[Flatten[Position[diagram,#]&/@hook]]];
Return[AllTrue[allhooks,balancedQ[#,inducedfilling[#]]&]];
];
variables=Apply[x,diagram,{1}];
flagconditions=Table[Apply[x,square]<=square[[1]],{square,diagram}];
columns=Table[Select[diagram,#[[2]]==p&],{p,DeleteDuplicates[(Reverse/@Sort[Reverse/@diagram])[[All,2]]]}];
strictcolumns=Flatten[Table[#[[1]]!=#[[2]]&/@Subsets[Apply[x,columns[[i]],{1}],{2}],{i,1,Length[columns]}]];
unbalancedfillings=variables/.Solve[Join[strictcolumns,flagconditions,Thread[variables>=1]],variables,Integers];
balancedlabelings=Select[unbalancedfillings,fillingbalancedQ];
matrices={};
Do[
A=ConstantArray[0,{n,n}];
Do[
A[[diagram[[i]][[1]],diagram[[i]][[2]]]]=bl[[i]];
,{i,1,Length[bl]}];
AppendTo[matrices,A];
,{bl,balancedlabelings}];
Return[matrices];
];

PartitionQ[lambda_]:=Module[{},
Return[VectorQ[lambda,IntegerQ[#]&&(Positive[#]||#==0)&]&&SortBy[lambda,-#&]==lambda];
];

PartitionToGrassmannianPermutation[lambda_]:=Module[{r,n,w},
r=Length[lambda];
w=Table[lambda[[r-i+1]]+i,{i,1,r}];
n=Max[w];
w=Join[w,Complement[Range[n],w]];
Return[w];
];

(*Options[SchurPolynomial]={Memoization->True};*)
SchurPolynomial[lambda_?PartitionQ,n_Integer?Positive,OptionsPattern[]]:=Module[{l,mu,delta},
l=Length[lambda];
If[n<l,Return[0]];
If[l<n,mu=Join[lambda,Table[0,{n-l}]],mu=lambda];
delta=Reverse[Range[0,n-1]];
Return[PolynomialReduce[Det[Table[x[i]^(mu[[j]]+delta[[j]]),{i,1,n},{j,1,n}]],{Det[Table[x[i]^delta[[j]],{i,1,n},{j,1,n}]]},Table[x[k],{k,1,n}]][[1,1]]];
];
SchurPolynomial[lambda_?PartitionQ]:=SchurPolynomial[lambda,Length[lambda]];

WordQ[list_?ListQ]:=Module[{},
Return[VectorQ[list,IntegerQ[#]&&#>=1&]];
];

WiringDiagram[inputword_?WordQ,n_Integer?Positive]:=Module[{word,points,strandpositions,colors,letter,higherstrand,lowerstrand,lines,leftlabels,permutation,rightlabels},
word=Reverse[inputword];
points=Reverse[Table[{{0,k}},{k,0,Max[word]}]];
strandpositions=Reverse[Table[{0,k},{k,0,Max[word]}]];
colors={Red,Orange,Yellow,Green,Blue,Cyan,Purple,Magenta,Brown,Pink,Gray,Black};
Do[
letter=word[[k]];
higherstrand=Position[strandpositions,{k-1,n-letter}][[1,1]];
lowerstrand=Position[strandpositions,{k-1,n-letter-1}][[1,1]];
strandpositions=strandpositions+Table[Switch[i,higherstrand,{1,-1},lowerstrand,{1,1},_,{1,0}],{i,1,n}];
Do[
AppendTo[points[[j]],strandpositions[[j]]];
,{j,1,n}]
,{k,1,Length[word]}];
lines=Flatten[Table[{Thick,colors[[Mod[i,Length[colors],1]]],Line[#]}&/@Partition[points[[i]],2,1],{i,1,Length[points]}],1];
leftlabels=Table[Style[Text[i,{-.5,n-i}],Large,FontSize->Large],{i,1,n}];
permutation=Reverse[Table[Position[strandpositions[[All,2]],l][[1,1]],{l,0,n-1}]];
rightlabels=Table[Style[Text[permutation[[i]],{Length[word]+.5,n-i}],Large,FontSize->Large],{i,1,n}];
Return[Graphics[Join[lines,leftlabels,rightlabels],PlotRange->{{-1,Length[word]+1},{-1/2,n-1/2}}]];
];
WiringDiagram[inputword_?WordQ]:=WiringDiagram[inputword,Max[inputword]+1];

Schubitope[w_?PermutationListQ]:=Module[{diagram,size,allSquares,word,phi,theta},
diagram=RotheDiagram[w];
size=Max[w];
allSquares=Tuples[Range[size],{2}];
word[c_,S_]:=DeleteCases[Table[If[!MemberQ[diagram,square]&&MemberQ[S,square[[1]]],"(",If[MemberQ[diagram,square]&&!MemberQ[S,square[[1]]],")",If[MemberQ[diagram,square]&&MemberQ[S,square[[1]]],"*"]]],{square,Select[allSquares,#[[2]]==c&]}],Null];
phi[c_,S_]:=Module[{count,wordcS,parentheses,occurences,stars},
count=0;
wordcS=word[c,S];
parentheses=Apply[StringJoin,DeleteCases[wordcS,"*"]];
occurences=StringCases[parentheses,"()"];
While[Length[occurences]>0,
count=count+Length[occurences];
parentheses=StringDelete[parentheses,"()"];
occurences=StringCases[parentheses,"()"];
];
stars=Length[StringCases[Apply[StringJoin,wordcS],"*"]];
Return[count+stars];
];
theta[S_]:=Sum[phi[c,S],{c,1,size}];
Return[theta];
];

SchubitopeDimension[w_?PermutationListQ]:=Module[{diagram,columns,columnindices,interval,intervals,intersectingpairs,joined},
diagram=RotheDiagram[w];
columns=Table[Select[diagram,#[[2]]==j&],{j,1,Length[w]}];
columnindices=columns[[All,All,1]];
interval[col_]:=If[col==Range[Length[col]],{},Range[Min[Complement[Range[Max[w]],col]],Max[col]]];
intervals=interval/@columnindices;
intervals=DeleteCases[intervals,{}];
intersectingpairs=Select[Tuples[Range[Length[intervals]],{2}],DeleteDuplicates[Sort[#]]==#&&IntersectingQ[intervals[[#[[1]]]],intervals[[#[[2]]]]]&];
While[intersectingpairs!={},joined=Union[intervals[[intersectingpairs[[1,1]]]],intervals[[intersectingpairs[[1,2]]]]];
intervals=Delete[intervals,{{intersectingpairs[[1,1]]},{intersectingpairs[[1,2]]}}];
intervals=Append[intervals,joined];
intersectingpairs=Select[Tuples[Range[Length[intervals]],{2}],DeleteDuplicates[Sort[#]]==#&&IntersectingQ[intervals[[#[[1]]]],intervals[[#[[2]]]]]&];
];
Return[Sum[Length[int]-1,{int,intervals}]];
];

ToZeroOneMatrix[diagram_?BoxListQ,{m_,n_}]:=Module[{A},
A=ConstantArray[0,{m,n}];
Do[A[[sq[[1]],sq[[2]]]]=1,{sq,diagram}];
Return[A];
];
ToZeroOneMatrix[diagram_?ListQ,{n_}]:=ToZeroOneMatrix[diagram,{n,n}];

DiagramMatrixQ[D_?MatrixQ]:=Module[{integer},
integer=MatrixQ[D,IntegerQ[#]&];
Return[integer];
];

KohnertMoveResults[diagram_?DiagramMatrixQ]:=Module[{m,n,squares,maxcolindices,rightmostsquares,newdiagrams,relevantcolumn,newrowindex,newdiagram},
m=Length[diagram];
n=Length[diagram[[1]]];
squares=Position[diagram,1];
maxcolindices=Max/@Table[Select[squares,#[[1]]==k&][[All,2]],{k,1,m}];
rightmostsquares=Select[Partition[Riffle[Range[m],maxcolindices],2],#[[2]]!=-\[Infinity]&];
newdiagrams={};
Do[
relevantcolumn=Select[squares,#[[2]]==rmsq[[2]]&&#[[1]]<=rmsq[[1]]&];
If[relevantcolumn[[All,1]]==Range[1,rmsq[[1]]],Continue[]];
newrowindex=Max[Complement[Range[1,rmsq[[1]]],relevantcolumn[[All,1]]]];
newdiagram=diagram;
newdiagram[[rmsq[[1]],rmsq[[2]]]]=0;
newdiagram[[newrowindex,rmsq[[2]]]]=1;
AppendTo[newdiagrams,newdiagram];
,{rmsq,rightmostsquares}];
Return[newdiagrams]
];

KohnertDiagrams[diagram_?DiagramMatrixQ]:=Module[{kohdiagrams,newdiagrams},
kohdiagrams={diagram};
newdiagrams=KohnertMoveResults[diagram];
While[newdiagrams!={},
kohdiagrams=DeleteDuplicates[Join[kohdiagrams,newdiagrams]];
newdiagrams=Flatten[KohnertMoveResults/@newdiagrams,1];
];
Return[kohdiagrams];
];

KohnertWeight[diagram_?DiagramMatrixQ]:=Module[{},
Return[Count[#,1]&/@diagram]
];

KohnertPolynomial[diagram_?DiagramMatrixQ]:=Module[{m,vars,weights},
m=Length[diagram];
vars=Array[x,m];
weights=KohnertWeight/@KohnertDiagrams[diagram];
Return[Total[Apply[Times,vars^#&/@weights,1]]];
];

RandomDiagram[m_Integer?Positive,n_Integer?Positive]:=RandomChoice[{0,1},{m,n}];

DrawKohnertDiagram[kdiagram_?DiagramMatrixQ]:=Module[{n,vertgridlines,horizgridlines,drawbox,drawcross,squares,crosses,drawing},
n=Length[kdiagram];
vertgridlines=Table[{Thick,Line[{{0,j},{n,j}}]},{j,0,n}];
horizgridlines=Table[{Thick,Line[{{j,0},{j,n}}]},{j,0,n}];
drawbox[{i_,j_}]:={Purple,Rectangle[#[[1]],#[[2]]]}&@({{#[[1]]-1,#[[2]]-1},{#[[1]],#[[2]]}}&@({#[[2]],n-#[[1]]+1}&@{i,j}));
drawcross[{i_,j_}]:={{Green,Thick,Line[{#[[1]],#[[2]]}]}&@({{#[[1]]-1,#[[2]]-1},{#[[1]],#[[2]]}}&@({#[[2]],n-#[[1]]+1}&@{i,j})),{Green,Thick,Line[{#[[1]]+{1,0},#[[2]]+{-1,0}}]}&@({{#[[1]]-1,#[[2]]-1},{#[[1]],#[[2]]}}&@({#[[2]],n-#[[1]]+1}&@{i,j}))};
squares=drawbox/@Join[Position[kdiagram,1],Position[kdiagram,2]];
crosses=Flatten[drawcross/@Position[kdiagram,2],1];
drawing=Graphics[Join[squares,crosses,horizgridlines,vertgridlines]]
];

GhostKohnertMoveResults[diagram_?DiagramMatrixQ]:=Module[{m,n,squares,deadsquares,maxcolindices,rightmostsquares,newdiagrams,relevantcolumn,zeros,exchangepos,newdiagram},
m=Length[diagram];
n=Length[diagram[[1]]];
squares=Position[diagram,1];
deadsquares=Position[diagram,2];
maxcolindices=Max/@Table[Select[Join[squares,deadsquares],#[[1]]==k&][[All,2]],{k,1,m}];
rightmostsquares=Select[Partition[Riffle[Range[m],maxcolindices],2],#[[2]]!=-\[Infinity]&];
rightmostsquares=Select[rightmostsquares,diagram[[#[[1]],#[[2]]]]==1&];
newdiagrams={};
Do[
relevantcolumn=Transpose[diagram][[rmsq[[2]],1;;rmsq[[1]]]];
If[!MemberQ[relevantcolumn,0],Continue[]];
zeros=Position[relevantcolumn,0];
exchangepos=Max[Flatten[zeros]];
If[MemberQ[relevantcolumn[[exchangepos;;rmsq[[1]]]],2],Continue[]];
newdiagram=diagram;
newdiagram[[rmsq[[1]],rmsq[[2]]]]=0;
newdiagram[[exchangepos,rmsq[[2]]]]=1;
AppendTo[newdiagrams,newdiagram];
newdiagram=diagram;
newdiagram[[rmsq[[1]],rmsq[[2]]]]=2;
newdiagram[[exchangepos,rmsq[[2]]]]=1;
AppendTo[newdiagrams,newdiagram];
,{rmsq,rightmostsquares}];
Return[DeleteDuplicates[newdiagrams]];
];

GhostKohnertDiagrams[diagram_?DiagramMatrixQ]:=Module[{kohdiagrams,newdiagrams},
kohdiagrams={diagram};
newdiagrams=GhostKohnertMoveResults[diagram];
While[newdiagrams!={},
kohdiagrams=DeleteDuplicates[Join[kohdiagrams,newdiagrams]];
newdiagrams=Flatten[GhostKohnertMoveResults/@newdiagrams,1];
];
Return[kohdiagrams];
];

GhostKohnertWeight[diagram_?DiagramMatrixQ]:=Module[{},
Return[(Count[#,1]&/@diagram)+(Count[#,2]&/@diagram)]
];

GhostKohnertPolynomial[diagram_?DiagramMatrixQ]:=Module[{m,vars,weights},
m=Length[diagram];
vars=Array[x,m];
weights=GhostKohnertWeight/@GhostKohnertDiagrams[diagram];
Return[Total[Apply[Times,vars^#&/@weights,1]]];
];

DrawDiagram[diagram_?BoxListQ,{m_,n_}]:=
Module[{vertgridlines,horizgridlines,squares,drawing},
horizgridlines=Table[{Thick,Line[{{0,j},{m,j}}]},{j,0,n}];
vertgridlines=Table[{Thick,Line[{{j,0},{j,n}}]},{j,0,m}];
squares={Purple,Rectangle[#[[1]],#[[2]]]}&/@({{#[[1]]-1,#[[2]]-1},{#[[1]],#[[2]]}}&/@({#[[2]],n-#[[1]]+1}&/@diagram));
drawing=Graphics[Join[squares,horizgridlines,vertgridlines]];
Return[drawing];
];
DrawDiagram[diagram_?BoxListQ,{n_}]:=DrawDiagram[diagram,{n,n}];
DrawDiagram[diagram_?BoxListQ]:=DrawDiagram[diagram,{Max[Flatten[diagram]]}];

ApplyKohnertMove[kdiagram_?DiagramMatrixQ,r_]:=Module[{rowrones,c,relevantcolumnslice,newrowindex,newdiagram},
rowrones=Flatten[Position[kdiagram[[r]],1]];
If[rowrones=={},Return[kdiagram]];
c=Max[rowrones];
relevantcolumnslice=Transpose[kdiagram][[c]][[1;;r]];
If[!MemberQ[relevantcolumnslice,0],Return[kdiagram]];
newrowindex=Max[Flatten[Position[relevantcolumnslice,0]]];
newdiagram=kdiagram;
newdiagram[[r,c]]=0;
newdiagram[[newrowindex,c]]=1;
Return[newdiagram];
];

ApplyKohnertSequence[kdiagram_?DiagramMatrixQ,seq_]:=Module[{diagram},
diagram=kdiagram;
Do[
diagram=ApplyKohnertMove[diagram,i];
,{i,seq}];
Return[diagram];
];

ValidKohnertMoveQ[diagram_?DiagramMatrixQ,r_]:=Module[{newdiagram},
newdiagram=ApplyKohnertMove[diagram,r];
Return[newdiagram!=diagram];
];

ValidKohnertSequenceQ[diagram_?DiagramMatrixQ,seq_]:=Module[{newdiagram,temp},
If[seq=={},Return[True]];
newdiagram=diagram;
Do[
temp=newdiagram;
newdiagram=ApplyKohnertMove[newdiagram,s];
If[temp==newdiagram,Return[False]];
,{s,seq}];
Return[True];
];

ValidKohnertMoves[kdiagram_?DiagramMatrixQ]:=Module[{},
Return[Select[Range[Length[kdiagram]],ValidKohnertMoveQ[kdiagram,#]&]];
];

ValidKohnertSequences[kdiagram_?DiagramMatrixQ]:=Module[{validmoves,sequences,move,nd,vks,newsequences},
validmoves=ValidKohnertMoves[kdiagram];
If[validmoves=={},Return[{{}}]];
sequences=Append[{#}&/@validmoves,{}];
Do[
move=validmoves[[i]];
nd=ApplyKohnertMove[kdiagram,move];
vks=ValidKohnertSequences[nd];
newsequences=Prepend[#,move]&/@vks;
sequences=Join[sequences,newsequences];
,{i,1,Length[validmoves]}];
Return[DeleteDuplicates[sequences]];
];

ValidRightmostLadderMoveQ[P_?DiagramMatrixQ,i_Integer?Positive]:=Module[{row,col,relevantcolumnslice,toprow,crosses1,crosses2},
row=P[[i]];
If[!MemberQ[row,1],Return[False]];
col=Max[Flatten[Position[row,1]]];
relevantcolumnslice=Transpose[P][[col]][[1;;i]];
If[!MemberQ[relevantcolumnslice,0],Return[False]];
toprow=Max[Flatten[Position[relevantcolumnslice,0]]];
If[P[[toprow,col+1]]==1,Return[False]];
crosses1=Transpose[P][[col]][[toprow+1;;i-1]];
crosses2=Transpose[P][[col+1]][[toprow+1;;i-1]];
If[crosses1!=Table[1,{Length[crosses1]}],Return[False]];
If[crosses2!=Table[1,{Length[crosses2]}],Return[False]];
Return[True];
];

ApplyRightmostLadderMove[P_?DiagramMatrixQ,i_Integer?Positive]:=Module[{row,col,relevantcolumnslice,toprow,crosses1,crosses2,newpipedream},
row=P[[i]];
If[!MemberQ[row,1],Return[P]];
col=Max[Flatten[Position[row,1]]];
relevantcolumnslice=Transpose[P][[col]][[1;;i]];
If[!MemberQ[relevantcolumnslice,0],Return[P]];
toprow=Max[Flatten[Position[relevantcolumnslice,0]]];
If[P[[toprow,col+1]]==1,Return[P]];
crosses1=Transpose[P][[col]][[toprow+1;;i-1]];
crosses2=Transpose[P][[col+1]][[toprow+1;;i-1]];
If[crosses1!=Table[1,{Length[crosses1]}],Return[P]];
If[crosses2!=Table[1,{Length[crosses2]}],Return[P]];
newpipedream=P;
newpipedream[[toprow,col+1]]=1;
newpipedream[[i,col]]=0;
Return[newpipedream];
];

ApplyRightmostLadderMoveSequence[P1_?DiagramMatrixQ,seq_?VectorQ]:=Module[{P},
P=P1;
Do[
P=ApplyRightmostLadderMove[P,i];
,{i,seq}];
Return[P];
];

ValidRightmostLadderMoves[P_?DiagramMatrixQ]:=Module[{},
Return[Select[Range[Length[P]],ValidRightmostLadderMoveQ[P,#]&]];
];

ValidRightmostLadderMoveSequenceQ[P1_?DiagramMatrixQ,seq_?VectorQ]:=Module[{P,flag},
If[seq=={},Return[True]];
P=P1;
flag=True;
Do[
If[!MemberQ[ValidRightmostLadderMoves[P],s],
flag=False;
Break[];
];
P=ApplyRightmostLadderMove[P,s];
,{s,seq}];
Return[flag];
];

ValidRightmostLadderMoveSequences[P_?DiagramMatrixQ]:=Module[{validmoves,sequences,move,npd,vrlms,newsequences},
validmoves=ValidRightmostLadderMoves[P];
If[validmoves=={},Return[{{}}]];
sequences=Append[{#}&/@validmoves,{}];
Do[
move=validmoves[[i]];
npd=ApplyRightmostLadderMove[P,move];
vrlms=ValidRightmostLadderMoveSequences[npd];
newsequences=Prepend[#,move]&/@vrlms;
sequences=Join[sequences,newsequences];
,{i,1,Length[validmoves]}];
Return[DeleteDuplicates[sequences]];
];

PipeDreamWeight[P_?PipeDreamQ]:=Module[{},
Return[Count[#,1]&/@P];
];

Options[Mitosis]={Reduced->True};
Mitosis[i_Integer?Positive,P_?PipeDreamQ,OptionsPattern[]]:=Module[{row,start,J,ones,newdiagrams,newdiagram,leftones,sums},
row=P[[i]];
start=Min[Flatten[Position[row,0]]];
J=Select[Range[start-1],P[[i+1,#]]!=1&];
ones=Flatten[Position[row,1]];
newdiagrams={};
Do[
newdiagram=P;
leftones=Select[ones,#<p&&MemberQ[J,#]&];
newdiagram[[i,p]]=0;
Do[
newdiagram[[i,k]]=0;
newdiagram[[i+1,k]]=1;
,{k,leftones}];
AppendTo[newdiagrams,newdiagram];
,{p,J}];
If[!OptionValue[Reduced],sums=Table[Sign[newdiagrams[[k]]+newdiagrams[[k+1]]],{k,1,Length[newdiagrams]-1}],sums={}];
Return[Join[newdiagrams,sums]];
];

FultonEssentialSet[w_]:=Module[{rothediagram},
rothediagram=RotheDiagram[w];
Return[Select[rothediagram,!MemberQ[rothediagram,#+{1,0}]&&!MemberQ[rothediagram,#+{0,1}]&]];
];

Compositions[n_Integer?Positive,{k_Integer?Positive}]:=DeleteDuplicates[Flatten[Permutations/@(PadLeft[#,k]&/@IntegerPartitions[n]),1]];
Compositions[n_Integer?Positive,k_Integer?Positive]:=DeleteDuplicates[Flatten[Compositions[n,{#}]&/@Range[k],1]];

FundamentalSlidePolynomial[alpha_?VectorQ]:=Module[{flat,dominatesQ,refinesQ,n,comps,sumover,variables},
If[alpha=={},Return[0]];
If[alpha==Table[0,{Length[alpha]}],Return[1]];
flat[a_]:=Select[a,#>0&];
dominatesQ[b_,a_]:=DeleteDuplicates[Table[Sum[b[[j]],{j,1,i}]>=Sum[a[[j]],{j,1,i}],{i,1,Length[a]}]]=={True};
refinesQ[b_,a_]:=Module[{psa,psb,positions},
psa=Accumulate[a];
psb=Accumulate[b];
If[!SubsetQ[psb,psa],Return[False]];
positions=Flatten[Position[psa,#]&/@psb];
Return[Sort[positions]==positions];
];
n=Length[alpha];
comps=Compositions[Total[alpha],{n}];
sumover=Select[comps,dominatesQ[#,alpha]&&refinesQ[flat[#],flat[alpha]]&];
variables=Table[x[i],{i,Range[n]}];
Return[Total[Table[Times@@sld,{sld,Table[variables^exp,{exp,sumover}]}]]];
];

FundamentalSlideExpansion[w_?PermutationListQ]:=Module[{sw,poly,comps,exp},
sw=SchubertPolynomial[w];
If[sw==1,Return[{{0}}]];
poly=sw;
comps={};
While[Exponents[poly,XVariables[w]]!={},
exp=Sort[Exponents[poly,XVariables[w]]][[1]];
poly=poly-FundamentalSlidePolynomial[exp];
AppendTo[comps,exp];
];
Return[comps];
];

MinimalPermutation[w_?PermutationListQ]:=Module[{max},
max=Max[Select[Range[Length[w]],w[[#]]!=#&]];
If[max==-\[Infinity],Return[{1}]];
Return[w[[1;;max]]];
];

Options[SchubertPolynomial]={Memoization->True};
SchubertPolynomial[w_?PermutationListQ,OptionsPattern[]]:=SchubertPolynomial[w,OptionValue[Memoization]];
SchubertPolynomial[w_,False]:=Module[{w0,u,word,g,i},w0=Reverse[Range[Length[w]]];
u=PermutationProduct[w0,InversePermutation[w]];
word=ReducedWord[u];
g=Apply[Times,Array[x[#]^(Length[w]-#)&,Length[w]-1]];
For[i=1,i<=Length[word],i++,g=DividedDifference[word[[i]],g];];
Return[Expand[g]];
];
SchubertPolynomial[w_,True]:=MemoizedSchubertPolynomial[MinimalPermutation[w]];

MemoizedSchubertPolynomial[w_]:=MemoizedSchubertPolynomial[w]=Module[{trans},
If[w==Range[Length[w]],Return[1]];
trans=TransitionRule[w];
Return[Expand[(x[trans[[1,1]]])MemoizedSchubertPolynomial[MinimalPermutation[trans[[1,2]]]]+Sum[MemoizedSchubertPolynomial[MinimalPermutation[u]],{u,trans[[2]]}]]];
];

Options[DoubleSchubertPolynomial]={Memoization->True};
DoubleSchubertPolynomial[w_?PermutationListQ,OptionsPattern[]]:=DoubleSchubertPolynomial[w,OptionValue[Memoization]];
DoubleSchubertPolynomial[w_,False]:=Module[{trans},
If[w==Range[Length[w]],Return[1]];
trans=TransitionRule[w];
Return[Expand[(x[trans[[1,1]]]-y[trans[[1,2]][[trans[[1,1]]]]])DoubleSchubertPolynomial[trans[[1,2]]]+Sum[DoubleSchubertPolynomial[u],{u,trans[[2]]}]]]
];
DoubleSchubertPolynomial[w_,True]:=MemoizedDoubleSchubertPolynomial[MinimalPermutation[w]];

MemoizedDoubleSchubertPolynomial[w_]:=MemoizedDoubleSchubertPolynomial[w]=Module[{trans},
If[w==Range[Length[w]],Return[1]];
trans=TransitionRule[w];
Return[Expand[(x[trans[[1,1]]]-y[trans[[1,2]][[trans[[1,1]]]]])MemoizedDoubleSchubertPolynomial[MinimalPermutation[trans[[1,2]]]]+Sum[MemoizedDoubleSchubertPolynomial[MinimalPermutation[u]],{u,trans[[2]]}]]];
];

Options[GrothendieckPolynomial]={Memoization->True};
GrothendieckPolynomial[w_?PermutationListQ,OptionsPattern[]]:=GrothendieckPolynomial[w,OptionValue[Memoization]];
GrothendieckPolynomial[w_,False]:=Module[{w0,u,word,g,i},
w0=Reverse[Range[Length[w]]];
u=PermutationProduct[w0,InversePermutation[w]];
word=ReducedWord[u];
g=Apply[Times,Array[x[#]^(Length[w]-#)&,Length[w]-1]];
For[i=1,i<=Length[word],i++,
g=IsobaricDividedDifference[word[[i]],g];
];
Return[Expand[g]];
];
GrothendieckPolynomial[w_,True]:=MemoizedGrothendieckPolynomial[MinimalPermutation[w]];

MemoizedGrothendieckPolynomial[w_]:=MemoizedGrothendieckPolynomial[w]=Module[{r,s,v,tterms,f,transition},
If[w==Sort[w],Return[1]];
r=Max[Descents[w]];
s=Max[Select[Range[r+1,Length[w]],w[[#]]<w[[r]]&]];
v=PermutationProduct[Transposition[r,s],w];
tterms=t[#,r]&/@Select[Range[r-1],NumberInversions[PermutationProduct[Transposition[#,r],v]]==NumberInversions[v]+1&];
f[set_]:=If[Length[set]==0,MemoizedGrothendieckPolynomial[MinimalPermutation[v]],(-1)^Length[set] MemoizedGrothendieckPolynomial[MinimalPermutation[PermutationProduct[PermutationProduct@@Reverse[(set/.{t[i_,j_]:>Transposition[i,j]})],v]]]];
transition=Total[f/@Subsets[tterms]];
Return[Expand[MemoizedGrothendieckPolynomial[MinimalPermutation[v]]+(x[r]-1)transition]];
];

Options[DoubleGrothendieckPolynomial]={Memoization->True};
DoubleGrothendieckPolynomial[w_?PermutationListQ,OptionsPattern[]]:=DoubleGrothendieckPolynomial[w,OptionValue[Memoization]];
DoubleGrothendieckPolynomial[w_?PermutationListQ,False]:=Module[{w0,u,word,g,i},
w0=Reverse[Range[Length[w]]];
u=PermutationProduct[w0,InversePermutation[w]];
word=ReducedWord[u];
g=Apply[Times,(x[#[[1]]]+y[#[[2]]]-x[#[[1]]]y[#[[2]]])&/@Select[Tuples[Range[Length[w]-1],{2}],#[[1]]+#[[2]]<=Length[w]&]];
For[i=1,i<=Length[word],i++,
g=IsobaricDividedDifference[word[[i]],g];
];
Return[Expand[g]];
];
DoubleGrothendieckPolynomial[w_,True]:=MemoizedDoubleGrothendieckPolynomial[MinimalPermutation[w]];

MemoizedDoubleGrothendieckPolynomial[w_]:=MemoizedDoubleGrothendieckPolynomial[w]=Module[{r,s,v,tterms,f,transition},
If[w==Sort[w],Return[1]];
r=Max[Descents[w]];
s=Max[Select[Range[r+1,Length[w]],w[[#]]<w[[r]]&]];
v=PermutationProduct[Transposition[r,s],w];
tterms=t[#,r]&/@Select[Range[r-1],NumberInversions[PermutationProduct[Transposition[#,r],v]]==NumberInversions[v]+1&];
f[set_]:=If[Length[set]==0,MemoizedDoubleGrothendieckPolynomial[MinimalPermutation[v]],(-1)^Length[set] MemoizedDoubleGrothendieckPolynomial[MinimalPermutation[PermutationProduct[PermutationProduct@@Reverse[(set/.{t[i_,j_]:>Transposition[i,j]})],v]]]];
transition=Total[f/@Subsets[tterms]];
Return[Expand[MemoizedDoubleGrothendieckPolynomial[MinimalPermutation[v]]+(1-x[r])(-1+y[v[[r]]])transition]];
];

(*
SchubertExpansion::failure="Output expansion equality check FAIL.";
Options[SchubertExpansion]={Memoization->True};
SchubertExpansion[poly_?IntegerQ,OptionsPattern[]]:=Return[{{Range[2]},{poly}}];
SchubertExpansion[poly_?IntegerQ,oldvars_?ListQ,OptionsPattern[]]:=Return[{Range[Max[1,Length[oldvars]]],{poly}}];
SchubertExpansion[poly_,OptionsPattern[]]:=SchubertExpansion[Expand[poly],Sort[Variables[Expand[poly]]],Memoization->OptionValue[Memoization]];
SchubertExpansion[oldpoly_,oldvars_?ListQ,OptionsPattern[]]:=Module[{poly,p,vars,exps,code,coeffs,codecoeff,perms,expcoeffs,w,sanitycheck},
poly=Expand[oldpoly];
If[poly==0,Return[{{},{}}]];
perms={};
expcoeffs={};
p=Expand[poly/.Table[oldvars[[i]]->x[i],{i,1,Length[oldvars]}]];
vars=Array[x,Length[vars]];
exps=Exponents[p,vars];
While[Length[exps]>0,
If[p==0,Break[]];
code=Sort[exps][[1]];
coeffs=Coefficients[p,vars];
codecoeff=coeffs[[Position[exps,code][[1,1]]]];
w=LehmerCodeToPermutation[code];
AppendTo[perms,w];
AppendTo[expcoeffs,codecoeff];
p=Expand[p-codecoeff*SchubertPolynomial[w,OptionValue[Memoization]]];
exps=Exponents[p,vars];
];
sanitycheck=Expand[(SchubertPolynomial/@perms).expcoeffs-poly];
If[NumberQ[sanitycheck]&&sanitycheck==0,
Return[{perms,expcoeffs}],Message[SchubertExpansion::failure];
];
];
*)

SchubertExpansion::failure="Output expansion equality check FAIL.";
Options[SchubertExpansion]={Memoization->True};
SchubertExpansion[poly_?IntegerQ,OptionsPattern[]]:=Return[{{Range[3]},{poly}}];
SchubertExpansion[poly_?IntegerQ,xvars_?ListQ,OptionsPattern[]]:=Return[{Range[Max[1,Length[xvars]]],{poly}}];
SchubertExpansion[poly_,OptionsPattern[]]:=Module[{vars,n},
vars=Variables[poly];
n=Max[Cases[vars,x[_]]/.x[a_]:>a];
Return[SchubertExpansion[Expand[poly],Array[x,n],Memoization->OptionValue[Memoization]]];
];
SchubertExpansion[oldpoly_,xvars_?ListQ,OptionsPattern[]]:=Module[{poly,p,exps,code,coeffs,codecoeff,perms,expcoeffs,w,n,paddedperms,sanitycheck},
poly=Expand[oldpoly];
If[poly==0, 
Return[{{},{}}];
];
perms={};
expcoeffs={};
p=poly;
exps=Exponents[p,xvars];
While[Length[exps]>0,
If[p==0,Break[]];
code=Sort[exps][[1]];
coeffs=Coefficients[p,xvars];
codecoeff=coeffs[[Position[exps,code][[1,1]]]];
w=LehmerCodeToPermutation[code];
AppendTo[perms,w];
AppendTo[expcoeffs,codecoeff];
p=Expand[p-codecoeff*SchubertPolynomial[w,OptionValue[Memoization]]];
exps=Exponents[p,xvars];
];
n=Max[Length/@perms];
paddedperms=Join[#,Complement[Range[n],#]]&/@perms;
sanitycheck=Expand[(SchubertPolynomial/@paddedperms).expcoeffs-poly];
If[NumberQ[sanitycheck]&&sanitycheck==0,
Return[{paddedperms,expcoeffs}],
Message[SchubertExpansion::failure]
];
];

KeyExpansion[w_?PermutationListQ]:=Module[{sw,poly,comps,exp},
sw=SchubertPolynomial[w];
If[sw==1,Return[{{0}}]];
poly=sw;
comps={};
While[Exponents[poly,XVariables[w]]!={},
exp=Sort[Exponents[poly,XVariables[w]]][[1]];
poly=poly-KeyPolynomial[exp];
AppendTo[comps,exp];
];
Return[comps];
];

NorthWestDiagramMatrixQ[matrix_]:=Module[{ones,pairs,goodpairs},
If[AnyTrue[Flatten[M],#!=0&&#!=1&],Return[False]];
ones=Position[matrix,1];
pairs=Tuples[ones,{2}];
goodpairs=Select[pairs,#[[1,1]]<#[[2,1]]&&#[[1,2]]>#[[2,2]]&];
Return[AllTrue[goodpairs,MemberQ[ones,{#[[1,1]],#[[2,2]]}]&]];
];

StronglySeparatedDiagramMatrixQ[M_]:=Module[{DD,pairs,TestPair},
If[AnyTrue[Flatten[M],#!=0&&#!=1&],Return[False]];
DD=Flatten[Position[#,1]]&/@Transpose[M];
pairs=Subsets[DD,{2}];
TestPair[pair_]:=Module[{int,comp1,comp2},
int=Intersection[pair[[1]],pair[[2]]];
comp1=Complement[pair[[1]],int];
comp2=Complement[pair[[2]],int];
If[comp1=={}||comp2=={},Return[True]];
Return[Max[comp1]<Min[comp2]||Max[comp2]<Min[comp1]];
];
Return[AllTrue[pairs,TestPair]];
];

Options[SchubertPolynomialQ]={Memoization->True};
SchubertPolynomialQ[oldpoly_,OptionsPattern[]]:=Module[{poly,vars,n,newvariables,f,exps,code,w,memoize,difference},
poly=Expand[oldpoly];
vars=Variables[poly];
n=Max[Cases[vars,x[_]]/.x[a_]:>a];
newvariables=Array[x,n];
f=poly;
exps=Exponents[f,Variables[f]];
code=Sort[exps][[1]];
w=LehmerCodeToPermutation[code];
memoize=OptionValue[Memoization];
difference=Expand[SchubertPolynomial[w,memoize]-f];
Return[Length[Exponents[difference,newvariables]]==0];
];

PercentageAvoidingDiagramMatrixQ[matrix_]:=Module[{ones,pairs,goodpairs},
ones=Position[matrix,1];
pairs=Tuples[ones,{2}];
goodpairs=Select[pairs,#[[1,1]]<#[[2,1]]&&#[[1,2]]>#[[2,2]]&];
Return[AllTrue[goodpairs,MemberQ[ones,{#[[1,1]],#[[2,2]]}]||MemberQ[ones,{#[[2,1]],#[[1,2]]}]&]];
];

SkewSchurPolynomial[lambda_?PartitionQ,mu1_?PartitionQ]:=Module[{mu,h,l,skewschur},
mu=PadRight[mu1,Length[lambda]];
h[k_,n_]:=Module[{vars,sols,polyvars},
If[k==0,Return[1]];
If[k<0,Return[0]];
vars=Array[t,n];
sols=vars/.Solve[Join[Thread[vars>=0],{Total[vars]==k}],vars,Integers];
If[sols==vars,sols={}];
polyvars=Array[x,n];
Return[Total[Times@@(polyvars^#)&/@sols]];
];
l=Length[lambda];
skewschur=Det[Table[h[lambda[[i]]-mu[[j]]-i+j,l],{i,1,l},{j,1,l}]];
Return[Expand[skewschur]];
];

MConvexSupportQ[f_]:=Module[{vars,exps},
If[IntegerQ[f],Return[True]];
vars=Variables[f];
exps=Exponents[f,vars];
Return[MConvexSetQ[exps]];
];

NormalizePolynomial[poly_]:=Module[{vars,monomials,exps,normalizers},
vars=Variables[poly];
monomials=MonomialList[poly];
exps=Flatten[Exponents[#,vars]&/@monomials,1];
normalizers=Apply[Times,Map[Factorial,exps,{2}],{1}];
Return[Sum[monomials[[i]]/normalizers[[i]],{i,1,Length[monomials]}]];
];

QuadraticFormToMatrix[f_,n_]:=Module[{indexvaluepairs,A},
If[IntegerQ[f]&&f==0,Return[ConstantArray[0,{n,n}]]];
indexvaluepairs=MonomialList[f]/.{a_ x[i_] x[j_]:>{i,j,a/2},x[i_] x[j_]:>{i,j,1/2},a_ x[i_]^2:>{i,i,a},x[i_]^2:>{i,i,1}};
A=ConstantArray[0,{n,n}];
Do[A[[ivp[[1]],ivp[[2]]]]=ivp[[3]];,{ivp,indexvaluepairs}];
Do[A[[s[[1]],s[[2]]]]=A[[s[[2]],s[[1]]]];,{s,Reverse/@Subsets[Range[n],{2}]}];
Return[A];
];

Options[LorentzianPolynomialQ]={CheckNonnegativity->True,CheckMConvexity->True};
LorentzianPolynomialQ[h_,OptionsPattern[]]:=Module[{oldvars,n,vars,f,degrees,homogeneous,coeffs,nonnegativecoeffs,d,g,A,mconvexsupport,increasingTuples,tuples,eigenvalues,spectral},
If[IntegerQ[h],Return[h>=0]];
oldvars=Variables[h];
n=Length[oldvars];
vars=Array[x,n];
f=h/.Thread[oldvars->vars];
degrees=Total/@Exponents[f,vars];
homogeneous=(Length[DeleteDuplicates[degrees]]==1);
If[!homogeneous,Return[False]];
coeffs=Coefficients[f,vars];
If[OptionValue[CheckNonnegativity],
nonnegativecoeffs=(DeleteDuplicates[NumberQ[#]&&#>=0&/@coeffs]=={True});
If[!nonnegativecoeffs,Return[False]];
];
If[OptionValue[CheckMConvexity],
mconvexsupport=MConvexSupportQ[f];
If[!mconvexsupport,Return[False]];
];
d=degrees[[1]];
If[d==1,Return[True]];
If[d==2,
g=f;
A=QuadraticFormToMatrix[g,n];
Return[Length[Select[Eigenvalues[A]//N,Positive]]<=1];
];
increasingTuples[l_List,m_]:=With[{iter={#2,#1,Length[l]}&@@@Partition[Prepend[Array[Unique[]&,m],1],2,1]},Part[l,#]&/@Flatten[Table@@{iter[[All,1]],Sequence@@iter},m-1]];
tuples=increasingTuples[Range[n],d-2];
spectral=True;
Do[
g=f;
Do[
g=D[g,x[i]];
,{i,Reverse[t]}];
A=QuadraticFormToMatrix[g,n];
eigenvalues=Eigenvalues[A]//N;
spectral=(Length[Select[eigenvalues,Positive[#]&]]<=1);
If[!spectral,Break[]];
,{t,tuples}];
Return[spectral];
];

HeckeReduce[inputword_]:=Module[{word,deletionpositions,deleterepeats,commutationpositions,commute,upbraidpositions,downbraidpositions,upbraid,downbraid,ubp,dbp,cmp,deletionsuccessful,newwords,words,stillgettingnewwords,upwords,downwords,commutewords,newerwords,max},
deletionpositions[w_]:=Select[Range[Length[w]-1],w[[#]]==w[[#+1]]&];
deleterepeats[w_]:=Module[{delposes,neww},
delposes=deletionpositions[w];
neww=w;
While[delposes!={},
neww=Delete[neww,delposes[[1]]];
delposes=deletionpositions[neww];
];
Return[neww];
];
commutationpositions[w_]:=Select[Range[Length[w]-1],Abs[w[[#]]-w[[#+1]]]>1&];
commute[w_,pos_]:=Join[w[[1;;pos-1]],{w[[pos+1]],w[[pos]]},w[[pos+2;;]]];
upbraidpositions[w_]:=Select[Range[Length[w]-2],w[[#]]==w[[#+1]]-1&&w[[#]]==w[[#+2]]&];
downbraidpositions[w_]:=Select[Range[Length[w]-2],w[[#]]==w[[#+1]]+1&&w[[#]]==w[[#+2]]&];
upbraid[w_,pos_]:=Join[w[[1;;pos-1]],{w[[pos]]+1,w[[pos+1]]-1,w[[pos+2]]+1},w[[pos+3;;]]];
downbraid[w_,pos_]:=Join[w[[1;;pos-1]],{w[[pos]]-1,w[[pos+1]]+1,w[[pos+2]]-1},w[[pos+3;;]]];
word=inputword;
word=deleterepeats[word];
Label[deletionsuccessful];
newwords={};
words={word};
Label[stillgettingnewwords];
Do[
ubp=upbraidpositions[w];
dbp=downbraidpositions[w];
cmp=commutationpositions[w];
upwords=upbraid[w,#]&/@ubp;
downwords=downbraid[w,#]&/@dbp;
commutewords=commute[w,#]&/@cmp;
newerwords=Select[DeleteDuplicates[Join[upwords,downwords,commutewords]],!MemberQ[words,#]&];
newwords=DeleteDuplicates[Join[newwords,newerwords]];
max=Max[Length[deletionpositions[#]]&/@newerwords];
If[max>0,Break[]];
,{w,words}];
If[max>0,
word=Select[newerwords,Length[deletionpositions[#]]>0&][[1]];
word=deleterepeats[word];
Goto[deletionsuccessful];
];
If[Length[newwords]>0,
words=Join[words,newwords];
newwords={};
Goto[stillgettingnewwords];
];
Return[words[[1]]];
];

PolynomialDegree[f_]:=Module[{exps,degree},
If[NumberQ[f]&&f==0,Return[-\[Infinity]]];
If[NumberQ[f],Return[0]];
exps=Exponents[f];
degree=Max[Total/@exps];
Return[degree];
];

RotheBPD[w_?PermutationListQ]:=Module[{P,j},
P=Table[0,{Length[w]},{Length[w]}];
For[j=1,j<=Length[w],j++,
P[[j,w[[j]]]]=1;
];
Return[P];
];

ASMQ[P_?MatrixQ]:=Module[{i,j, sum, Valid},
Valid=True;
For[i=1,i<=Length[P],i++,
For[j=1,j<=Length[P],j++,
Valid=If[P[[i,j]]!=0 && P[[i,j]]!=1&&P[[i,j]]!=-1, False, Valid];
];
];
For[i=1,i<=Length[P],i++,
sum=0;
For[j=1,j<=Length[P],j++,
sum=sum+P[[i,j]];
Valid=If[sum!=0 && sum!=1, False, Valid];
];
];
For[j=1,j<=Length[P],j++,
sum=0;
For[i=1,i<=Length[P],i++,
sum=sum+P[[i,j]];
Valid=If[sum!=0 && sum!=1, False, Valid];
];
];
Return[Valid];
];

DroopMove[P_?MatrixQ, i0_, j0_, imax_, jmax_]:=Module[{i,j,newP},
newP=P;
newP[[i0,j0]]=0;
newP[[i0,jmax]]=1;
newP[[imax,j0]]=1;
newP[[imax,jmax]]=-1;
Return[newP];
];

DroopMoveQ[P_?MatrixQ, i0_, j0_, imax_, jmax_]:=Module[{i,j, Valid},
Valid=True;
Valid=If[P[[i0,j0]]!=1,False, Valid];
For[i=i0,i<=imax,i++,
For[j=j0,j<=jmax,j++,
Valid=If[P[[i,j]]!=0 && (i!=i0||j!=j0), False, Valid];
];
];
Valid=If[ASMQ[DroopMove[P,i0,j0,imax,jmax]],Valid,False];
Return[Valid];
];

ReducedBumplessPipeDreams[w_?PermutationListQ]:=Module[{P0, P, BPD, NewBPD, i, j, imax, jmax,  k, n},
P0=RotheBPD[w];
BPD={P0};
NewBPD=BPD;
n=0;
While[n!=Length[BPD],
n=Length[BPD];
For[j=1,j<=Length[w]-1,j++,
For[i=1,i<=Length[w]-1,i++,
For[jmax=j+1,jmax<=Length[w],jmax++,
For[imax=i+1,imax<=Length[w],imax++,
For[k=1,k<=Length[BPD],k++,
P=If[DroopMoveQ[BPD[[k]],i,j,imax,jmax], DroopMove[BPD[[k]],i,j,imax,jmax],P0];
If[!MemberQ[NewBPD,P]&&ASMQ[P],AppendTo[NewBPD,P]];
DeleteDuplicates[NewBPD];
];
BPD=NewBPD;
];
];
];
];
];
Return[BPD];
];

BPDtoPermutation[P_?ASMQ]:=Module[{i,j, k,l, a, b, word, Perm, DownElbowPositions, UpElbowPositions,VerticalPositions,HorizontalPositions, CrossPositions, temp},
DownElbowPositions=Position[P,1];
UpElbowPositions=Position[P,-1];
VerticalPositions={};
HorizontalPositions={};
For[k=1,k<=Length[P],k++,
l=1;
While[P[[Length[P]+1-l,k]]==0,
If[l==1&&k==1,VerticalPositions={{Length[P]+1-l,k}},VerticalPositions=Append[VerticalPositions,{Length[P]+1-l,k}]];
l++;
];
];
For[k=1,k<=Length[UpElbowPositions],k++,
l=1;
While[P[[UpElbowPositions[[k,1]]-l,UpElbowPositions[[k,2]]]]==0,
VerticalPositions=Append[VerticalPositions,{UpElbowPositions[[k,1]]-l,UpElbowPositions[[k,2]]}];
l++;
];
];
For[k=1,k<=Length[DownElbowPositions],k++,
l=1;
While[DownElbowPositions[[k,2]]+l<=Length[P]&&P[[DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l]]==0,
If[l==1&&k==1,HorizontalPositions={{DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l}},HorizontalPositions=Append[HorizontalPositions,{DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l}]];
l++;
];
];
CrossPositions=Intersection[HorizontalPositions,VerticalPositions];
word={};
For[j=1,j<=Length[P],j++,
For[i=Length[P],i>=1,i--,
If[MemberQ[CrossPositions,{i,j}],
word=Append[word,Sum[P[[a,b]],{a,1,i},{b,1,j}]-1]
];
];
];
Perm=Range[Length[P]];
For[k=1,k<=Length[word],k++,
temp=Perm[[word[[k]]+1]];
Perm[[word[[k]]+1]]=Perm[[word[[k]]]];
Perm[[word[[k]]]]=temp;
];
Return[Perm];
];

Options[DrawBumplessPipeDream]={PipeColor->Green,Labels->True};
DrawBumplessPipeDream[P_?ASMQ,OptionsPattern[]]:=Module[{vertgridlines,horizgridlines,DownElbowPositions, UpElbowPositions,VerticalPositions,HorizontalPositions,hline,hlines,vline,vlines,downelbow,downelbows,upelbow,upelbows,Pcrosspositions,changeofcoordinates,Crossing,crosses,lowerlabels,perm,rightlabels,picture,k,l},
vertgridlines=Table[{Thick,Line[{{0,j},{Length[P],j}}]},{j,0,Length[P]}];
horizgridlines=Table[{Thick,Line[{{j,0},{j,Length[P]}}]},{j,0,Length[P]}];
DownElbowPositions=Position[P,1];
UpElbowPositions=Position[P,-1];
VerticalPositions={};
HorizontalPositions={};
For[k=1,k<=Length[P],k++,
l=1;
While[P[[Length[P]+1-l,k]]==0,
If[l==1&&k==1,VerticalPositions={{Length[P]+1-l,k}},VerticalPositions=Append[VerticalPositions,{Length[P]+1-l,k}]];
l++;
];
];
For[k=1,k<=Length[UpElbowPositions],k++,
l=1;
While[P[[UpElbowPositions[[k,1]]-l,UpElbowPositions[[k,2]]]]==0,
VerticalPositions=Append[VerticalPositions,{UpElbowPositions[[k,1]]-l,UpElbowPositions[[k,2]]}];
l++;
];
];
For[k=1,k<=Length[DownElbowPositions],k++,
l=1;
While[DownElbowPositions[[k,2]]+l<=Length[P]&&P[[DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l]]==0,
If[l==1&&k==1,HorizontalPositions={{DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l}},HorizontalPositions=Append[HorizontalPositions,{DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l}]];
l++;
];
];
changeofcoordinates[{i_,j_}]:={j-1/2,Length[P]+1/2-i};
If[HorizontalPositions!={},
hline[{i_,j_}]:={{Thick,OptionValue[PipeColor],Line[{{i-1/2,j},{i+1/2,j}}]}};
hlines=hline/@(changeofcoordinates/@HorizontalPositions);
,hlines={};
];
If[VerticalPositions!={},
vline[{i_,j_}]:={{Thick,OptionValue[PipeColor],Line[{{i,j-1/2},{i,j+1/2}}]}};
vlines=vline/@(changeofcoordinates/@VerticalPositions);
,vlines={};
];
If[DownElbowPositions!={},
downelbow[{i_,j_}]:={{Thick,OptionValue[PipeColor],Circle[{i+1/2,j-1/2},1/2,{Pi/2,Pi}]}};
downelbows=downelbow/@(changeofcoordinates/@DownElbowPositions);
,downelbows={};
];
If[UpElbowPositions!={},
upelbow[{i_,j_}]:={{Thick,OptionValue[PipeColor],Circle[{i-1/2,j+1/2},1/2,{0,-Pi/2}]}};
upelbows=upelbow/@(changeofcoordinates/@UpElbowPositions);
,upelbows={};
];
lowerlabels=Table[{Text[Style[i,Large],{i-1/2,-1/5}]},{i,1,Length[P]}];
perm=BPDtoDemazureProduct[P];
rightlabels=Table[{Text[Style[perm[[i]],Large],{Length[P]+1/5,Length[P]+1/2-i}]},{i,1,Length[P]}];
If[OptionValue[Labels],
picture=Graphics[Join[vertgridlines,horizgridlines,vlines, hlines, upelbows,downelbows,lowerlabels,rightlabels]],
picture=Graphics[Join[vertgridlines,horizgridlines,vlines, hlines, upelbows,downelbows]]
];
Return[picture];
];

kDroop1[P_?MatrixQ, i0_, j0_, imax_, jmax_]:=Module[{i,j,elbow, newP},
elbow=1;
For[j=j0,j<=jmax-1,j++,
elbow=If[P[[imax,j]]==1, j, elbow];
];
newP=P;
newP[[i0,j0]]=0;
newP[[imax,j0]]=1;
newP[[imax,elbow]]=0;
newP[[i0,elbow]]=1;
Return[newP];
];

KDroop1Q[P_?MatrixQ, i0_, j0_, imax_, jmax_]:=Module[{i,j, Valid, counter},
Valid=True;
Valid=If[P[[i0,j0]]!=1,False, Valid];
Valid=If[P[[imax,jmax]]!=-1,False, Valid];
For[i=i0,i<=imax-1,i++,
For[j=j0,j<=jmax,j++,
Valid=If[P[[i,j]]!=0 && (i!=i0||j!=j0), False, Valid];
];
];
counter=0;
For[j=j0,j<=jmax-1,j++,
counter=If[P[[imax,j]]==1, counter+1, counter];
Valid=If[P[[imax,j]]==-1 , False, Valid];
];
Valid=If[counter==1, Valid, False];
Valid=If[ASMQ[kDroop1[P,i0,j0,imax,jmax]],Valid,False];
Return[Valid];
];

kDroop2[P_?MatrixQ, i0_, j0_, imax_, jmax_]:=Module[{i,j,elbow, newP},
elbow=1;
For[i=i0,i<=imax-1,i++,
elbow=If[P[[i,jmax]]==1, i, elbow];
];
newP=P;
newP[[i0,j0]]=0;
newP[[i0,jmax]]=1;
newP[[elbow,jmax]]=0;
newP[[elbow,j0]]=1;
Return[newP];
];

KDroop2Q[P_?MatrixQ, i0_, j0_, imax_, jmax_]:=Module[{i,j, Valid, counter},
Valid=True;
Valid=If[P[[i0,j0]]!=1,False, Valid];
Valid=If[P[[imax,jmax]]!=-1,False, Valid];
For[i=i0,i<=imax,i++,
For[j=j0,j<=jmax-1,j++,
Valid=If[P[[i,j]]!=0 && (i!=i0||j!=j0), False, Valid];
];
];
counter=0;
For[i=i0,i<=imax-1,i++,
counter=If[P[[i,jmax]]==1, counter+1, counter];
Valid=If[P[[i,jmax]]==-1 , False, Valid];
];
Valid=If[counter==1, Valid, False];
Valid=If[ASMQ[kDroop2[P,i0,j0,imax,jmax]],Valid,False];
Return[Valid];
];

BumplessPipeDreams[w_?PermutationListQ]:=Module[{P0, P, BPD, NewBPD, i, j, imax, jmax, k, n},
P0=RotheBPD[w];
BPD={P0};
NewBPD=BPD;
n=0;
While[n!=Length[BPD],
n=Length[BPD];
For[j=1,j<=Length[w]-1,j++,
For[i=1,i<=Length[w]-1,i++,
For[jmax=j+1,jmax<=Length[w],jmax++,
For[imax=i+1,imax<=Length[w],imax++,
For[k=1,k<=Length[BPD],k++,
P=If[DroopMoveQ[BPD[[k]],i,j,imax,jmax], DroopMove[BPD[[k]],i,j,imax,jmax],P0];
If[!MemberQ[NewBPD,P]&&ASMQ[P],AppendTo[NewBPD,P]];
DeleteDuplicates[NewBPD];
];
BPD=NewBPD;
];
];
NewBPD=BPD;
For[jmax=j+1,jmax<=Length[w],jmax++,
For[imax=i+1,imax<=Length[w],imax++,
For[k=1,k<=Length[BPD],k++,
P=If[KDroop1Q[BPD[[k]],i,j,imax,jmax], kDroop1[BPD[[k]],i,j,imax,jmax],P0];
If[!MemberQ[NewBPD,P]&&ASMQ[P],AppendTo[NewBPD,P]];
DeleteDuplicates[NewBPD];
];
BPD=NewBPD;
];
];
NewBPD=BPD;
For[jmax=j+1,jmax<=Length[w],jmax++,
For[imax=i+1,imax<=Length[w],imax++,
For[k=1,k<=Length[BPD],k++,
P=If[KDroop2Q[BPD[[k]],i,j,imax,jmax], kDroop2[BPD[[k]],i,j,imax,jmax],P0];
If[!MemberQ[NewBPD,P]&&ASMQ[P],AppendTo[NewBPD,P]];
DeleteDuplicates[NewBPD];
];
BPD=NewBPD;
];
];

];
];
];
Return[BPD];
];

PivotQ[P_?ASMQ, a_, b_, i_, j_]:=Module[{k,l, output},
output=True;
output=If[P[[i,j]]!= 1, False, output];
output=If[i<=a &&j<=b, output, False];
For[k=i,k<=a,k++,
For[l=j,l<=b,l++,
output=If[P[[k,l]]!= 0&&{k,l}!={i,j}&&{k,l}!={a,b}, False, output];
];
];
Return[output];
];

RemovableNeg1Q[P_?ASMQ, a_, b_]:=Module[{k,l, n, output, pivots},
output=If[P[[a,b]]!= -1, False, True];
pivots={};
For[k=1,k<=a,k++,
For[l=b,l>=1,l--,
If[PivotQ[P,a,b,k,l]&&{k,l}!={a,b},pivots=Append[pivots,{k,l}]
];
];
];
For[n=1,n<=Length[pivots]-1,n++,
For[k=pivots[[n,1]],k<=a,k++,
For[l=pivots[[n+1,2]],l<=b,l++,
output=If[P[[k,l]]== -1&&{k,l}!={a,b}, False, output];
];
];
];

Return[output];
];

BPDtoDemazureProduct[P_?ASMQ]:=Module[{i,j, k,l, a, b, n,  Perm, temp, pivots},
temp=P;
For[a=1, a<=Length[temp]-1, a++,
For[b=1, b<=Length[temp]-1, b++,
If[RemovableNeg1Q[temp,a,b], 
pivots={};
For[k=1,k<=a,k++,
For[l=b,l>=1,l--,
If[PivotQ[temp,a,b,k,l]&&{k,l}!={a,b},
pivots=Append[pivots,{k,l}];
];
];
];
For[n=1,n<=Length[pivots]-1,n++,
For[k=pivots[[n,1]],k<=a,k++,
For[l=pivots[[n+1,2]],l<=b,l++,
temp[[k,l]]=0;
];
];
];
For[n=1,n<=Length[pivots]-1,n++,
temp[[pivots[[n,1]],pivots[[n+1,2]]]]=1;
];

];
];
];
Perm=BPDtoPermutation[temp];
Return[Perm];
];

MarkedBumplessPipeDreams[w_?PermutationListQ]:=Module[{P0, P, BPD, NewBPD, i, j, imax, jmax,  k, mBPD, U, shaded},
BPD=BumplessPipeDreams[w];
mBPD={};
For[i=1, i<=Length[BPD], i++,
U=Position[BPD[[i]],-1];
shaded=Subsets[U];
For[j=1, j<=Length[shaded], j++,
AppendTo[mBPD, {BPD[[i]],shaded[[j]]}];
];
];
Return[mBPD];
];

Options[DrawMarkedBumplessPipeDream]={PipeColor->Green,Labels->True};
DrawMarkedBumplessPipeDream[{P_?ASMQ, s_},OptionsPattern[]]:=Module[{vertgridlines,horizgridlines,DownElbowPositions, UpElbowPositions,VerticalPositions,HorizontalPositions,hline,hlines,vline,vlines,downelbow,downelbows,upelbow,upelbows,Pcrosspositions,changeofcoordinates,Crossing,crosses,lowerlabels,perm,rightlabels,picture,k,l, shade, shades},
vertgridlines=Table[{Thick,Line[{{0,j},{Length[P],j}}]},{j,0,Length[P]}];
horizgridlines=Table[{Thick,Line[{{j,0},{j,Length[P]}}]},{j,0,Length[P]}];
DownElbowPositions=Position[P,1];
UpElbowPositions=Position[P,-1];
VerticalPositions={};
HorizontalPositions={};
For[k=1,k<=Length[P],k++,
l=1;
While[P[[Length[P]+1-l,k]]==0,
If[l==1&&k==1,VerticalPositions={{Length[P]+1-l,k}},VerticalPositions=Append[VerticalPositions,{Length[P]+1-l,k}]];
l++;
];
];
For[k=1,k<=Length[UpElbowPositions],k++,
l=1;
While[P[[UpElbowPositions[[k,1]]-l,UpElbowPositions[[k,2]]]]==0,
VerticalPositions=Append[VerticalPositions,{UpElbowPositions[[k,1]]-l,UpElbowPositions[[k,2]]}];
l++;
];
];
For[k=1,k<=Length[DownElbowPositions],k++,
l=1;
While[DownElbowPositions[[k,2]]+l<=Length[P]&&P[[DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l]]==0,
If[l==1&&k==1,HorizontalPositions={{DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l}},HorizontalPositions=Append[HorizontalPositions,{DownElbowPositions[[k,1]],DownElbowPositions[[k,2]]+l}]];
l++;
];
];
changeofcoordinates[{i_,j_}]:={j-1/2,Length[P]+1/2-i};
If[HorizontalPositions!={},
hline[{i_,j_}]:={{Thick,OptionValue[PipeColor],Line[{{i-1/2,j},{i+1/2,j}}]}};
hlines=hline/@(changeofcoordinates/@HorizontalPositions);
,hlines={};
];
If[VerticalPositions!={},
vline[{i_,j_}]:={{Thick,OptionValue[PipeColor],Line[{{i,j-1/2},{i,j+1/2}}]}};
vlines=vline/@(changeofcoordinates/@VerticalPositions);
,vlines={};
];
If[DownElbowPositions!={},
downelbow[{i_,j_}]:={{Thick,OptionValue[PipeColor],Circle[{i+1/2,j-1/2},1/2,{Pi/2,Pi}]}};
downelbows=downelbow/@(changeofcoordinates/@DownElbowPositions);
,downelbows={};
];
If[UpElbowPositions!={},
upelbow[{i_,j_}]:={{Thick,OptionValue[PipeColor],Circle[{i-1/2,j+1/2},1/2,{0,-Pi/2}]}};
upelbows=upelbow/@(changeofcoordinates/@UpElbowPositions);
,upelbows={};
];
If[s!={},
shade[{i_,j_}]:={LightGray,Rectangle[{i-1/2,j-1/2},{i+1/2,j+1/2}]};
shades=shade/@(changeofcoordinates/@s);
];
lowerlabels=Table[{Text[Style[i,Large],{i-1/2,-1/5}]},{i,1,Length[P]}];
perm=BPDtoDemazureProduct[P];
rightlabels=Table[{Text[Style[perm[[i]],Large],{Length[P]+1/5,Length[P]+1/2-i}]},{i,1,Length[P]}];
If[s=={},
If[OptionValue[Labels],
picture=Graphics[Join[vertgridlines,horizgridlines, vlines, hlines, upelbows,downelbows,lowerlabels,rightlabels]],
picture=Graphics[Join[vertgridlines,horizgridlines, vlines, hlines, upelbows,downelbows]]
];
,
If[OptionValue[Labels],
picture=Graphics[Join[shades, vertgridlines,horizgridlines, vlines, hlines, upelbows,downelbows,lowerlabels,rightlabels]],
picture=Graphics[Join[shades, vertgridlines,horizgridlines, vlines, hlines, upelbows,downelbows]]
];
];
Return[picture];
];

End[];

SetAttributes[{
Memoization,
Permutations,
PermutationListForm,
PipeDreamQ,
PipeDreamToPermutation,
PipeDreams,
DrawPipeDream,
Inversions,
NumberInversions,
ReducedPipeDreamQ,
ReducedPipeDreams,
GrothendieckPolynomial,
BetaGrothendieckPolynomial,
DoubleGrothendieckPolynomial,
BetaDoubleGrothendieckPolynomial,
SchubertPolynomial,
Descents,
Ascents,
Transposition,
TransitionRule,
DoubleSchubertPolynomial,
Letter,
ReducedWords,
ReducedWord,
DividedDifference,
IsobaricDividedDifference,
LehmerCode,
LehmerCodeToPermutation,
OneFixedDominantQ,
CoreRegion,
RotheDiagram,
BottomReducedPipeDream,
TopReducedPipeDream,
AvoidsPattern,
GeneralizedPermutahedronZVector,
GeneralizedPermutahedronInequalities,
GeneralizedPermutahedronYVector,
MConvexSetQ,
Coefficients,
Exponents,
DominantQ,
GrassmannianQ,
VexillaryQ,
DrawRotheDiagram,
DemazureDifference,
DemazureLascouxDifference,
KeyPolynomial,
HomogeneousComponents,
SymmetricPolynomialQ,
RankMatrix,
BruhatOrderLessEqualQ,
MultiplicityPolynomials,
XVariables,
XBVariables,
XYVariables,
XYBVariables,
BruhatOrderLessEqualCoverQ,
SkylineDiagram,
BalancedLabelings,
PartitionQ,
SchurPolynomial,
WiringDiagram,
Schubitope,
SchubitopeDimension,
ToZeroOneMatrix,
KohnertMoveResults,
KohnertDiagrams,
KohnertWeight,
KohnertPolynomial,
RandomDiagram,
DrawKohnertDiagram,
GhostKohnertMoveResults,
GhostKohnertDiagrams,
GhostKohnertWeight,
GhostKohnertPolynomial,
DrawDiagram,
ApplyKohnertMove,
ApplyKohnertSequence,
ValidKohnertMoveQ,
ValidKohnertSequenceQ,
ValidKohnertMoves,
ValidKohnertSequences,
ValidRightmostLadderMoveQ,
ApplyRightmostLadderMove,
ApplyRightmostLadderMoveSequence,
ValidRightmostLadderMoves,
ValidRightmostLadderMoveSequenceQ,
ValidRightmostLadderMoveSequences,
PipeDreamWeight,
Mitosis,
FultonEssentialSet,
FundamentalSlidePolynomial,
Compositions,
FundamentalSlideExpansion,
MinimalPermutation,
KeyExpansion,
SchubertExpansion,
NorthWestDiagramMatrixQ,
StronglySeparatedDiagramMatrixQ,
SchubertPolynomialQ,
PercentageAvoidingDiagramMatrixQ,
SkewSchurPolynomial,
MConvexSupportQ,
NormalizePolynomial,
QuadraticFormToMatrix,
LorentzianPolynomialQ,
CheckNonnegativity,
CheckMConvexity,
HeckeReduce,
PolynomialDegree,
RotheBPD,
ASMQ,
DroopMove,
DroopMoveQ,
ReducedBumplessPipeDreams,
BPDtoPermutation,
DrawBumplessPipeDream,
PipeColor,
kDroop1,
KDroop1Q,
kDroop2,
KDroop2Q,
BumplessPipeDreams,
PivotQ,
RemovableNeg1Q,
BPDtoDemazureProduct,
MarkedBumplessPipeDreams,
DrawMarkedBumplessPipeDream
},{Protected,ReadProtected}
];

EndPackage[];
