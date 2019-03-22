#!/bin/bash

changelog=""
function gen_changelog()
{
    local firstline=$( grep -n "^chameleon-" ChangeLog | head -n 1 | cut -d ':' -f 1 )
    firstline=$(( firstline + 2 ))
    #echo $firstline
    local lastline=$( grep -n "^chameleon-" ChangeLog | head -n 2 | tail -n 1 | cut -d ':' -f 1 )
    lastline=$(( lastline - 1 ))
    #echo $lastline

    changelog="Changes:\n"
    for i in `seq $firstline $lastline`
    do
        local line=$( head -n $i ChangeLog | tail -n 1 )
        changelog="$changelog$line\\n"
        #echo $line
    done

    changelog="${changelog}\n__WARNING__: Download the source archive by clicking on the link __Download release__ above, please do not consider the link Source code to get all submodules.\n"
}

release=""
function get_release()
{
    local firstline=$( grep -n "^chameleon-" ChangeLog | head -n 1 | cut -d ':' -f 1 )
    release=$( head -n $firstline ChangeLog | tail -n 1 | sed 's/chameleon\-//' )
}

# Get the release name through the branch name, and through the ChangeLog file.
# Both have to match to be correct
RELEASE_NAME=`echo $CI_COMMIT_REF_NAME | cut -d - -f 2`
get_release

if [ -z "$RELEASE_NAME" -o -z "$release" -o "$RELEASE_NAME" != "$release" ]
then
    echo "Commit name $RELEASE_NAME is different from ChangeLog name $release"
    exit 1
fi

# extract the change log from ChangeLog
gen_changelog
echo $changelog

# generate the archive
wget https://raw.githubusercontent.com/Kentzo/git-archive-all/master/git_archive_all.py
mv git_archive_all.py git-archive-all
chmod +x git-archive-all
./git-archive-all --force-submodules chameleon-$RELEASE_NAME.tar.gz

# upload the source archive
GETURL=`echo curl --request POST --header \"PRIVATE-TOKEN: $RELEASE_TOKEN\" --form \"file=\@chameleon-$RELEASE_NAME.tar.gz\" https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/uploads`
MYURL=`eval $GETURL | jq .url | sed "s#\"##g"`

# Try to remove the release if it already exists
curl --request DELETE --header "PRIVATE-TOKEN: $RELEASE_TOKEN" https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/releases/v$RELEASE_NAME

# create the release and the associated tag
COMMAND=`echo curl --header \"Content-Type: application/json\" --header \"PRIVATE-TOKEN: $RELEASE_TOKEN\" \
  --data \'{ \"name\": \"v$RELEASE_NAME\", \
            \"tag_name\": \"v$RELEASE_NAME\", \
            \"ref\": \"$CI_COMMIT_REF_NAME\", \
            \"description\": \"$changelog\", \
            \"assets\": { \"links\": [{ \"name\": \"Download release\", \"url\": \"$CI_PROJECT_URL$MYURL\" }] } }\' \
  --request POST https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/releases`
eval $COMMAND
