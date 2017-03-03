for i in *.fits.gz
#for i in 0589p000.fits 3584p121.fits
	do
		keep=0
		python wcs_wise.py $i
		polyid ../boss_survey.ply $i.coords $i.out
		sed -i '1d' $i.out
		while IFS=" " read -r val1 val2 val3
			do
				if [[ $val3 ]] 
					then keep=1
				fi
			done < "$i.out"

		if (($keep == 0))
			then mv $i outboss/.

		elif (($keep == 1))
			then mv $i inboss/. 
		fi
		rm -f $i.out $i.coords
		echo $i
	done
	
