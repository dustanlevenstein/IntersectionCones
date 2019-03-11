
LoadPackage("hecke");

SizeScreen([ 80, 43 ]); # calibrates computers.

True := true;
False := false;

settings := ReadAsFunction("settings.txt")();

# Generates an order function which sorts by p-cores first, so blocks are organized.
generate_order_function := function(p)
	local order_function;
	order_function := function(mu, nu)
		local ecore_mu, ecore_nu;
		ecore_mu := ECore(p, mu);
		ecore_nu := ECore(p, nu);
		if ecore_mu < ecore_nu then
			return true;
		elif ecore_mu = ecore_nu then
			return LengthLexicographic(mu, nu);
		else
			return false;
		fi;
	end;;
	return order_function;
end;;
ordering := generate_order_function(2); # Just for example. This should be changed when p is selected.

debug_print := function(p)
	Print(p);
	end;;

nonreduced_dual_cone_gens_raw := function(p, n, P_p, higher_order_cones)
	local k, l_p, proj, e, P_e, P_etrans, PB_e, QB_e, l_e, vec, temp_rec, p_part, i, entry;
	# First compute P, B, Q, B' for e=p.
	# P_p := MatrixDecompositionMatrix(DecompositionMatrix(Specht(p), n));
	k := Size(P_p);
	l_p := Size(P_p[1]);

	p_part := IdentityMat(l_p); # The positive part of the dual cone.

	e := p*p;
	for P_e in higher_order_cones do
		# P_e := MatrixDecompositionMatrix(DecompositionMatrix(Specht(e), n));
		l_e := Size(P_e[1]);
		P_etrans := TransposedMatMutable(P_e);
		temp_rec := BaseSteinitzVectors(IdentityMat(k), P_etrans); # Just an arbitrary way to extend the basis.
		Append(P_etrans, temp_rec.factorspace);
		PB_e := TransposedMatImmutable(P_etrans);
		QB_e := Inverse(PB_e);
		i := 1;
		while i <= k do
			if i <= l_e then
				Add(p_part, QB_e[i]*P_p);
			else
				# Verify pi annihilates this row (just a sanity check).
				for entry in QB_e[i]*P_p do
					if entry <> 0 then
						Print("Fuck.\n");
					fi;
				od;
			fi;
			i := i+1;
		od;
		e := e*p;
	od;
	Sort(p_part);
	i := 2;
	while i <= Length(p_part) do
		if p_part[i] = p_part[i-1] then
			Remove(p_part, i);
		else
			i := i+1;
		fi;
	od;
	return rec(P_p := P_p, p_part := p_part);
end;;

nonreduced_dual_cone_gens_params := function(p, n)
	local P_p, P_es, e;

	# First compute P, B, Q, B' for e=p.
	P_p := MatrixDecompositionMatrix(DecompositionMatrix(Specht(p), n, ordering));
	P_es := [];
	e := p*p;
	while e <= n do
		Add(P_es, MatrixDecompositionMatrix(DecompositionMatrix(Specht(e), n, ordering)));
		e := e*p;
	od;
	return rec(P_p := P_p, P_es := P_es);
end;;


nonreduced_dual_cone_gens := function(p, n)
	local r;
	r := nonreduced_dual_cone_gens_params(p, n);
	return nonreduced_dual_cone_gens_raw(p, n, r.P_p, r.P_es);
end;;

LoadPackage( "ctbllib" );


# Use stbllib instead of Hecke. Then I process the blocks in python because that's easier.
dump_true_cones := function(p, n)
	local sn, t, charparams, decmat, filename;
	ordering := generate_order_function(p);
	sn := Concatenation("S", ViewString(n));
	t := CharacterTable(sn);
	charparams := List(CharacterParameters(t), x -> x[2]);
	decmat := List(DecompositionMatrix(t mod p));
	SortParallel(charparams, decmat, ordering);
	filename := JoinStringsWithSeparator(
		["../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/true_decomposition_matrices/p_", 
		ViewString(p),
		"n_",
		ViewString(n),
		".txt"], "");
	PrintTo(filename, "return ");
	AppendTo(filename, rec( parts := charparams, decmat := decmat));
end;;

generate_partitions := function(p, n)
	local i, partitions, p_cores, previous, partitions_preg, partitions_p2reg, partitions_l3, partitions_l3_preg, partitions_l3_p2reg,
		partitions_eregs, partitions_l3_eregs, index, e;

	partitions := Reversed(ERegularPartitions(Specht(0), n));
	p_cores := [];
	previous := 0;
	Sort(partitions, ordering);
	for i in partitions do
		i := ECore(p, i);
		if previous = 0 then
			Add(p_cores, i);
			previous := i;
		elif previous <> i then
			Add(p_cores, i);
			previous := i;
		fi;
	od;
	partitions_preg := Reversed(ERegularPartitions(Specht(p), n));
	partitions_p2reg := Reversed(ERegularPartitions(Specht(p^2), n));
	partitions_eregs := [];
	e := p^2;
	i := 1;
	while e <= n do
		Add(partitions_eregs, Reversed(ERegularPartitions(Specht(e), n)));
		Sort(partitions_eregs[i], ordering);
		e := e*p;
		i := i+1;
	od;
	Sort(partitions_preg, ordering);
	Sort(partitions_p2reg, ordering);
	partitions_l3 := [];
	partitions_l3_preg := [];
	partitions_l3_p2reg := [];
	partitions_l3_eregs := [];
	for i in partitions_eregs do
		Add(partitions_l3_eregs, []);
	od;
	for i in partitions do
		if Length(i) <= 3 then
			Add(partitions_l3, i);
			if Position(partitions_preg, i) <> fail then
				Add(partitions_l3_preg, i);
				if Position(partitions_p2reg, i) <> fail then
					Add(partitions_l3_p2reg, i);
					index := 1;
					while index <= Length(partitions_eregs) do
						if Position(partitions_eregs[index], i) <> fail then
							Add(partitions_l3_eregs[index], i);
						fi;
						index := index+1;
					od;
				fi;
			fi;
		fi;
	od;
	return rec(partitions := partitions, partitions_preg := partitions_preg, partitions_p2reg := partitions_p2reg,
		partitions_l3 := partitions_l3, partitions_l3_preg := partitions_l3_preg, partitions_l3_p2reg := partitions_l3_p2reg,
		p_cores := p_cores, partitions_eregs := partitions_eregs, partitions_l3_eregs := partitions_l3_eregs);
end;;



dump_raw_cones := function(p, n)
	local parts, r, filename;
	ordering := generate_order_function(p);
	
	parts := generate_partitions(p, n);
	r := nonreduced_dual_cone_gens_params(p, n);
	filename := JoinStringsWithSeparator(
		["../../Dropbox/UCLA/Research/current_quarter/IntersectionConeUtils/actually_raw_cones/p_", 
			ViewString(p),
			"n_",
			ViewString(n),
			".txt"], "");
	PrintTo(filename, "return ");
	AppendTo(filename, rec( P_p := r.P_p, P_es := r.P_es, parts := parts));
end;;

dump_all_raw_cones_given_p := function(p)
	local n;
	n := 1;
	while true do
		debug_print(ViewString(n));
		debug_print("\n");
		dump_raw_cones(p, n);
		n := n+1;
	od;
end;;

dump_all_raw_cones_given_p := function(p, n)
	while true do
		debug_print(ViewString(n));
		debug_print("\n");
		dump_raw_cones(p, n);
		n := n+1;
	od;
end;;
