#find . -type f | xargs sed -i 's/and Freie Universitaet Berlin 2002-2022/and Freie Universitaet Berlin 2002-2022/g'
## For all files and on MacOS: 
find . -type f ! -path "./.git/*" -exec grep -q "Freie Universitaet Berlin 2002-20" {} \; -exec sed -i '' -e 's/and Freie Universitaet Berlin 2002-2022/and Freie Universitaet Berlin 2002-2022/g' {} \;
