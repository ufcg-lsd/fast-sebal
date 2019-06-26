#pragma once

#include "types.h"
#include "pixel_reader.h"

/**
 * @brief  Configures a TIFF based on a second TIFF.
 * @param  new_tif: TIFF to be configured.
 * @param  base_tif: TIFF used to provide the configurations.
 */
void setup(TIFF* new_tif, TIFF* base_tif);

/**
 * @brief  Verifies if a TIFF was open correctly. 
 * @param  tif: TIFF to be verified
 * @throws Throw an error with exit code 1 if the TIFF isn't open.
 */
void check_open_tiff(TIFF* tif);

/**
 * @brief  Reads the values of a line in a TIFF saving them into an array.
 * @param  tif: TIFF who line should be read.
 * @param  tif_line[]: Array where the data will be saved.
 * @param  line: Number of the line to be read.
 * @throws Throw an error with exit code 3 if the read couldn't be done.
 */
void read_line_tiff(TIFF* tif, double tif_line[], int line);

/**
 * @brief  Reads the values of a line in a TIFF saving them into an array.
 * @param  tif: TIFF who line should be read.
 * @param  tif_line: image data ref
 * @param  line: Number of the line to be read.
 * @throws Throw an error with exit code 3 if the read couldn't be done.
 */
void read_line_tiff(TIFF* tif, tdata_t tif_line, int line);

/**
 * @brief  Reads the value contained in a specific position of a TIFF.
 * @param  tif: TIFF who value should be read.
 * @param  col: Number of the column to be read.
 * @param  line: Number of the line to be read.
 * @throws Throw an error with exit code 3 if the read couldn't be done. 
 */
double read_position_tiff(TIFF* tif, int col, int line);

/**
 * @brief  Writes values from an array to a specific line in a TIFF.
 * @param  tif: TIFF who line should be written.
 * @param  tif_line[]: Array containing the values to be written.
 * @param  line: Number of the line to be read.
 * @throws Throw an error with exit code 4 if the write couldn't be done.
 */
void write_line_tiff(TIFF* tif, double tif_line[], int line);

/**
 * @brief  Closes open TIFFs.
 * @param  tiffs[]: Array containing opened tiffs to be closed.
 * @param  quant_tiffs: Length of the array or number of tiffs.
 */
void close_tiffs(TIFF* tiffs[], int quant_tiffs);

/*
The following definitions are from The art of computer programming by Knuth
*/

/**
 * @brief  Determines if a and b are approximately equals based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if they are approximately equals, and FALSE otherwise.
 */
bool approximatelyEqual(double a, double b);

/**
 * @brief  Determines if a and b are essentially equals based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if they are essentially equals, and FALSE otherwise.
 */
bool essentiallyEqual(double a, double b);

/**
 * @brief  Determines if a is definitely greater than b based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if a is definitely greater than b, and FALSE otherwise.
 */
bool definitelyGreaterThan(double a, double b);

/**
 * @brief  Determines if a is definitely less than b based on a epsilon.
 * @param  a: First value.
 * @param  b: Second value.
 * @retval TRUE if a is definitely less than b, and FALSE otherwise.
 */
bool definitelyLessThan(double a, double b);