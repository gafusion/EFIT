#!/bin/shell
/f/python/hp/bin/gmake clean
/f/python/hp/bin/gmake EFIT
echo "efitd6565d built"
/f/python/hp/bin/gmake clean
/f/python/hp/bin/gmake 'grid=129' EFIT
echo "efitd129d built"
/f/python/hp/bin/gmake clean
/f/python/hp/bin/gmake 'grid=257' EFIT
echo "efitd257d built"
/f/python/hp/bin/gmake clean
/f/python/hp/bin/gmake 'grid=513' EFIT
echo "efitd513d built"
