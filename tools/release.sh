#!/usr/bin/env bash
###
#
#  @file release.sh
#  @copyright 2013-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Script to generate the release when pushing a branch and tag of the same name
#
#  @version 1.2.0
#  @author Florent Pruvost
#  @author Mathieu Faverge
#  @date 2022-02-22
#
###

#
# Steps to do the release:
#    - Update information in the code (see update_release.sh)
#    - Update the ChangeLog
#    - Push the hash on solverstack as:
#          - a tag named vx.x.x
#          - a branch named release-x.x.x (will trigger the CI to generate the release)
#
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

    changelog="$changelog\nWARNING: Download the source archive by clicking on the link __Download release__ above, please do not consider the automatic Source code links as they are missing the submodules.\n"
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

# generate the archive
wget https://raw.githubusercontent.com/Kentzo/git-archive-all/master/git_archive_all.py
python3 git_archive_all.py --force-submodules chameleon-$RELEASE_NAME.tar.gz

# upload the source archive
GETURL=`echo curl --request POST --header \"PRIVATE-TOKEN: $RELEASE_TOKEN\" --form \"file=\@chameleon-$RELEASE_NAME.tar.gz\" https://gitlab.inria.fr/api/v4/projects/$CI_PROJECT_ID/uploads`
MYURL=`eval $GETURL | jq .url | sed "s#\"##g"`

# extract the change log from ChangeLog
gen_changelog
echo $changelog

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
