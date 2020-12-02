thisdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
thisfile=$(basename ${BASH_SOURCE[0]})
echo ${thisdir}
echo ${thisfile}

echo ${thisdir}/${thisfile}
