idxs=('17924y' '17927b' '17933b' '17934u' '17949r' '17960r' '17965r')
echo $idxs
for idx in "${idxs[@]}"; do
	echo "$idx"
	#mkdir ./loc_fiber_deci_EDT/${idx}
	mkdir ./loc_fiber_deci_transfer_EDT/${idx}
	echo "LOCFILES ./obs_deci_transfer/${idx}_picks.hyp HYPO71 ../DAS_location/time_punto_medio/layer ./loc_fiber_deci_transfer_EDT/${idx}/PhaseNet" >> run_punto_medio/NLLoc-PS.in
        NLLoc run_punto_medio/ISNet_nll.run	
done
