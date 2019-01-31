#!/bin/bash
# Calculate Hi-C map resolution
# calculation were followed https://github.com/aidenlab/juicer/edit/master/misc/calculate_map_resolution.sh

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [map files]

Description
	sort.map file has following column
	1. id
	2. chr 1
	3. position 1
	4. direction 1
	5. map quality 1
	6. fragment id 1
	7. fragment location 1
	8. repeat 1
	9. chr 2
	10. position 2
	11. direction 2
	12. map quality 2
	13. fragment id 2
	14. fragment location 2
	15. repeat 2
	list the all map files and make space between different files.

	-h, --help
		show help

	-v, --version
		show version

	-L, --length
		total chromosome length (default : 3095693983 for hg19. Specify value for other length)
		length information could be obtain from http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics

	-c, --count [count file]
		read count file for each 50bp. It will creat from map file. If file is not empty, file information will use resolution estimation. If empty, read count will calculate
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvL:c:
LONG=help,version,length:,count:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? -ne 0 ]]; then
	exit 2
fi
eval set -- "$PARSED"

while true; do
	case "$1" in
		-h|--help)
			get_usage
			exit 1
			;;
		-v|--version)
			get_version
			exit 1
			;;
		-L|--length)
			total="$2"
			shift 2
			;;
		-c|--count)
			coveragename="$2"
			shift 2
			;;			
		--)
			shift
			break
			;;
		*)
			echo "Programming error"
			exit 3
			;;
	esac
done

DIR_LIB=$(dirname $0)
TIME_STAMP=$(date +"%Y-%m-%d")

[ ! -n "${coveragename}" ] && echo "Please specify count file name" && exit 1
total=${total:-3095693983}

MAP_FILES=$@

# Create 50bp coverage vector
if [ ! -s $coveragename ]
then
    awk '{
      if ($5>0&&$12>0&&$6!=$13)
        {
        chr1=0;
        chr2=0;

        chr1=$2; 
        chr2=$9;
        if (chr1!=0&&chr2!=0)
        {
         val[chr1 " " int($3/50)*50]++
         val[chr2 " " int($10/50)*50]++
        }
      }
   }
   END{
     for (i in val)
     {
       print i, val[i]
     }
   }' "$MAP_FILES" > $coveragename
fi


# threshold is 80% of total bins
binstotal=$(( $total / 50 ))
threshold=$(( $binstotal * 4 ))
threshold=$(( $threshold / 5 ))

echo -ne "."
newbin=50
bins1000=$(awk '$3>=1000{sum++}END{if (sum == 0) print 0; else print sum}' $coveragename)
lowrange=$newbin


# find reasonable range with big jumps
while [ $bins1000 -lt $threshold ]
do
    lowrange=$newbin
    newbin=$(( $newbin + 1000 ))
    echo -ne "."
    bins1000=$(awk -v x=$newbin '{ 
      val[$1 " " int($2/x)*x]=val[$1 " " int($2/x)*x]+$3
    }
    END { 
      for (i in val) { 
        if (val[i] >= 1000) {
          count++
        } 
     } 
     print count
   }' $coveragename )
    binstotal=$(( $total / $newbin ))
    threshold=$(( $binstotal * 4 ))
    threshold=$(( $threshold / 5 ))
done

# at this point, lowrange failed but newbin succeeded
# thus the map resolution is somewhere between (lowrange, newbin]
midpoint=$(( $newbin - $lowrange ))
midpoint=$(( $midpoint / 2 ))
midpoint=$(( $midpoint + $lowrange ))
# now make sure it's a factor of 50 (ceil)
midpoint=$(( $midpoint + 49 ))
midpoint=$(( $midpoint / 50 ))
midpoint=$(( $midpoint * 50 ))

# binary search
while [ $midpoint -lt $newbin ]
do
    # echo -ne "."
    bins1000=$(awk -v x=$midpoint '{ 
      val[$1 " " int($2/x)*x]=val[$1 " " int($2/x)*x]+$3
    }
    END { 
      for (i in val) { 
        if (val[i] >= 1000) {
          count++
        } 
     } 
     print count
   }' $coveragename )
    binstotal=$(( $total / $midpoint ))
    threshold=$(( $binstotal * 4 ))
    threshold=$(( $threshold / 5 ))
    if [ $bins1000 -lt $threshold ]
    then
	lowrange=$midpoint;
	# at this point, lowrange failed but newbin succeeded
	midpoint=$(( $newbin - $lowrange ))
	midpoint=$(( $midpoint / 2 ))
	midpoint=$(( $midpoint + $lowrange ))
	# now make sure it's a factor of 50 (ceil)
	midpoint=$(( $midpoint + 49 ))
	midpoint=$(( $midpoint / 50 ))
	midpoint=$(( $midpoint * 50 ))
    else
	newbin=$midpoint;
	# at this point, lowrange failed but newbin succeeded
	midpoint=$(( $newbin - $lowrange ))
	midpoint=$(( $midpoint / 2 ))
	midpoint=$(( $midpoint + $lowrange ))
	# now make sure it's a factor of 50
	midpoint=$(( $midpoint + 49 ))
	midpoint=$(( $midpoint / 50 ))
	midpoint=$(( $midpoint * 50 ))
    fi
done

echo -e "\nThe map resolution is $newbin"