# Run

- Landsat 5 and 7
./run b1 b2 b3 b4 b5 b6 b7 bqa mtl raster_elevation output_path

- Landsat 8
./run b2 b3 b4 b5 b6 b7 b10 bqa mtl raster_elevation output_path

# Exit code

- 1 is open tiff problem
- 2 is open file problem (metadata or station data)
- 3 is read tiff problem
- 4 is write tiff problem
- 5 is shadow problem
- 6 is sensor problem
- 7 is type tiff problem (pixel reader)
