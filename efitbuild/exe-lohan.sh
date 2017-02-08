make clean
make EFIT
mv efitd6565dpgf90 efitd6565d
echo "efitd6565d built"
make clean
make 'grid=129' EFIT
mv efitd129dpgf90 efitd129d
echo "efitd129d built"
make clean
make 'grid=257' EFIT
mv efitd257dpgf90 efitd257d
echo "efitd257d built"
make clean
make 'grid=513' EFIT
mv efitd513dpgf90 efitd513d
echo "efitd513d built"
