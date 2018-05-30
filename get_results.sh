root_dir=$(pwd)

for folder in $root_dir/*
do 
	file_dir="$folder/result"
	for file in $file_dir/*
	do
		echo "$file"
		sed -n -e 30p $file
	done
	echo " "
done