#find . -type f | xargs sed -i 's/and Freie Universitaet Berlin 2002-2020/and Freie Universitaet Berlin 2002-2020/g'
## For all files and on MacOS: 
LC_ALL=C find . -type f ! -path "./.git/*" -exec grep -q "Freie Universitaet Berlin 2002-20" {} \; -exec sed -i '' -e 's/and Freie Universitaet Berlin 2002-2020/and Freie Universitaet Berlin 2002-2020/g' {} \;
