find . -type f | xargs sed -i 's/and Freie Universitaet Berlin 2002-20../and Freie Universitaet Berlin 2002-2018/g'
## For all files and on MacOS: find . -type f ! -path "./.git/*" -exec sed -i '' -e 's/and Freie Universitaet Berlin 2002-20../and Freie Universitaet Berlin 2002-2018/g' {} \;
