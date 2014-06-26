errc=0;

for m in $( seq 4 128 ); do
	echo -n ${m}

	hm=$(( m / 2 ))
	for b in $( seq 2 ${hm} ); do

		for t in $( seq 2 ${b} ); do

			for i in $( seq 1 10 ); do
				../bin/rlchol ${m} ${b} ${t} 1 2 &> /dev/null
				if [ $? -ne 0 ]; then
					echo -n"error:  "
					echo "${m} ${b} ${t} 1 2" 
					let errc=errc+1;
					exit 1
				fi
			done
			echo -n "."
		done
	done
	echo 
done

echo "found ${errc} errors"
