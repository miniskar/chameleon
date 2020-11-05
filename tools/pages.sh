#!/usr/bin/env bash

CHAMELEON_SRC_DIR=${CHAMELEON_SRC_DIR:-$PWD}

mkdir tmp_fig
cd tmp_fig

## need to generate figures from last benchmarks
# get the most recent chameleon commit date for which benchs have been performed
commit_sha=`curl -X GET "https://elasticsearch.bordeaux.inria.fr/hiepacs-chameleon_perf/_search?pretty" -H 'Content-Type: application/json' -d'
{
  "query": { "match_all": {} },
  "sort": [
    { "Commit_date_chameleon": "desc" }
  ],
  "size": 1
}
' | grep "Commit_sha_chameleon" | awk '{ print $3" "$4}' |sed -e "s#,##g" |sed -e "s#\"##g" |sed -e "s# ##g"`
echo $commit_sha

# generate the csv file from elasticsearch for the given chameleon commit
python3 ${CHAMELEON_SRC_DIR}/tools/bench/jube/get_result.py -e https://elasticsearch.bordeaux.inria.fr -t hiepacs -p chameleon -c $commit_sha

# generate the figures
Rscript ${CHAMELEON_SRC_DIR}/tools/bench/jube/GenFigures.R

# add the performance chapter to the doc. we need the performances files to add
# this chapter that are not commited in the sources so that this chapter is not
# here by default
cat >> ${CHAMELEON_SRC_DIR}/doc/orgmode/users_guide.org.in <<EOF
* Chameleon Performances on PlaFRIM
Chameleon commit: *$commit_sha*.
#+INCLUDE: @CMAKE_CURRENT_SOURCE_DIR@/chapters/performances.org
EOF

cd ..

## Build the doc
VERSION=${VERSION:-pages}
mkdir -p build-$VERSION
cd build-$VERSION

cmake $CHAMELEON_SRC_DIR -DCHAMELEON_ENABLE_DOC=ON
make doc -j5

cd ..
mv build-$VERSION/doc/orgmode          public/
mv build-$VERSION/doc/doxygen/out/html public/doxygen
mv tmp_fig/*.png public/

cd public/
if [ -f users_guide.html ]
then
    ln -sfn users_guide.html index.html
else
    echo -e "ERROR: missing users_guide.html file"
    exit 1
fi

# Change the width of the page in the CSS file
sed -i -e "s#max-width:800px#max-width:1800px#" org-html-themes/styles/readtheorg/css/readtheorg.css
