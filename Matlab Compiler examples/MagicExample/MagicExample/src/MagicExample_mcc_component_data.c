/*
 * MATLAB Compiler: 4.6 (R2007a)
 * Date: Wed Aug 08 14:46:54 2007
 * Arguments: "-B" "macro_default" "-o" "MagicExample" "-W" "main" "-d"
 * "C:\tamar\Matlab\MagicExample\MagicExample\src" "-T" "link:exe" "-v"
 * "C:\tamar\Matlab\MagicExample\magicsquare.m" 
 */

#include "mclmcr.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_MagicExample_session_key[] = {
        '3', '7', 'E', '6', 'F', '9', 'A', 'E', 'F', '0', '0', '2', 'B', '9',
        '5', '9', '1', 'D', '4', '4', '1', '2', 'F', '5', '8', 'A', '7', '9',
        '4', '4', 'B', 'E', '1', 'B', '7', '2', '9', '6', '3', 'D', '3', 'D',
        '2', 'C', 'E', 'B', '7', '9', '0', '4', 'B', '3', '7', '1', 'B', 'F',
        '8', '4', '7', '9', '1', 'F', '1', '1', '6', '3', '1', '4', '7', 'C',
        '1', '4', 'A', 'D', 'F', '4', '1', '6', '4', 'D', '0', 'C', '8', 'D',
        'C', '9', 'C', '8', '1', 'B', '4', '0', '0', '6', '9', '6', '9', 'B',
        '1', 'C', 'C', '9', 'D', '9', 'B', '3', '7', '7', 'D', 'B', '9', '6',
        '0', '9', '0', '9', '2', 'D', '4', '2', '8', '1', '9', 'D', 'D', '2',
        'B', 'B', 'B', 'A', '2', 'F', 'B', 'E', '6', '6', '5', 'B', '1', '5',
        'B', '7', 'E', '9', 'A', '8', '1', 'C', '1', '5', 'D', '4', 'F', '1',
        '1', '7', '5', '9', '1', '6', '4', '9', '6', '1', 'D', '9', '1', '0',
        '4', '7', '8', '4', 'C', 'B', '3', '6', 'C', '6', '4', '9', '1', 'C',
        'E', '1', 'A', '9', '0', 'E', '8', '2', '6', 'C', '9', '4', '2', '8',
        '0', 'E', '1', '3', '7', '8', '8', 'A', 'E', '0', '3', '4', '7', '0',
        '2', '7', '8', 'F', 'B', '5', 'D', 'F', 'D', 'E', 'E', '4', '0', '9',
        '1', '0', 'E', 'E', 'F', '8', '9', '2', 'A', 'E', 'B', 'F', '8', '6',
        '6', '3', '0', '0', '8', '6', 'B', 'F', '9', 'C', '9', 'D', '9', '3',
        'C', '0', '6', '9', '\0'};

const unsigned char __MCC_MagicExample_public_key[] = {
        '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9',
        '2', 'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1',
        '0', '1', '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B',
        '0', '0', '3', '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1',
        '0', '0', 'C', '4', '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3',
        'A', '5', '2', '0', '6', '5', '8', 'F', '6', 'F', '8', 'E', '0', '1',
        '3', '8', 'C', '4', '3', '1', '5', 'B', '4', '3', '1', '5', '2', '7',
        '7', 'E', 'D', '3', 'F', '7', 'D', 'A', 'E', '5', '3', '0', '9', '9',
        'D', 'B', '0', '8', 'E', 'E', '5', '8', '9', 'F', '8', '0', '4', 'D',
        '4', 'B', '9', '8', '1', '3', '2', '6', 'A', '5', '2', 'C', 'C', 'E',
        '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4', 'D', '0', '8', '5',
        'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2', 'E', 'D', 'E',
        '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6', '3', '7',
        '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E', '6',
        '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
        '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1',
        'B', 'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9',
        '9', '0', '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0',
        'B', '6', '1', 'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B',
        '5', '8', 'F', 'C', '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6',
        'E', 'B', '7', 'E', 'C', 'D', '3', '1', '7', '8', 'B', '5', '6', 'A',
        'B', '0', 'F', 'A', '0', '6', 'D', 'D', '6', '4', '9', '6', '7', 'C',
        'B', '1', '4', '9', 'E', '5', '0', '2', '0', '1', '1', '1', '\0'};

static const char * MCC_MagicExample_matlabpath_data[] = 
    { "MagicExample/", "toolbox/compiler/deploy/", "tamar/Matlab/",
      "tamar/Matlab/functions toolbox/", "tamar/Matlab/draw toolbox/",
      "$TOOLBOXMATLABDIR/general/", "$TOOLBOXMATLABDIR/ops/",
      "$TOOLBOXMATLABDIR/lang/", "$TOOLBOXMATLABDIR/elmat/",
      "$TOOLBOXMATLABDIR/elfun/", "$TOOLBOXMATLABDIR/specfun/",
      "$TOOLBOXMATLABDIR/matfun/", "$TOOLBOXMATLABDIR/datafun/",
      "$TOOLBOXMATLABDIR/polyfun/", "$TOOLBOXMATLABDIR/funfun/",
      "$TOOLBOXMATLABDIR/sparfun/", "$TOOLBOXMATLABDIR/scribe/",
      "$TOOLBOXMATLABDIR/graph2d/", "$TOOLBOXMATLABDIR/graph3d/",
      "$TOOLBOXMATLABDIR/specgraph/", "$TOOLBOXMATLABDIR/graphics/",
      "$TOOLBOXMATLABDIR/uitools/", "$TOOLBOXMATLABDIR/strfun/",
      "$TOOLBOXMATLABDIR/imagesci/", "$TOOLBOXMATLABDIR/iofun/",
      "$TOOLBOXMATLABDIR/audiovideo/", "$TOOLBOXMATLABDIR/timefun/",
      "$TOOLBOXMATLABDIR/datatypes/", "$TOOLBOXMATLABDIR/verctrl/",
      "$TOOLBOXMATLABDIR/codetools/", "$TOOLBOXMATLABDIR/helptools/",
      "$TOOLBOXMATLABDIR/winfun/", "$TOOLBOXMATLABDIR/demos/",
      "$TOOLBOXMATLABDIR/timeseries/", "$TOOLBOXMATLABDIR/hds/",
      "$TOOLBOXMATLABDIR/guide/", "$TOOLBOXMATLABDIR/plottools/",
      "toolbox/local/", "toolbox/compiler/" };

static const char * MCC_MagicExample_classpath_data[] = 
    { "" };

static const char * MCC_MagicExample_libpath_data[] = 
    { "" };

static const char * MCC_MagicExample_app_opts_data[] = 
    { "" };

static const char * MCC_MagicExample_run_opts_data[] = 
    { "" };

static const char * MCC_MagicExample_warning_state_data[] = 
    { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_MagicExample_component_data = { 

    /* Public key data */
    __MCC_MagicExample_public_key,

    /* Component name */
    "MagicExample",

    /* Component Root */
    "",

    /* Application key data */
    __MCC_MagicExample_session_key,

    /* Component's MATLAB Path */
    MCC_MagicExample_matlabpath_data,

    /* Number of directories in the MATLAB Path */
    39,

    /* Component's Java class path */
    MCC_MagicExample_classpath_data,
    /* Number of directories in the Java class path */
    0,

    /* Component's load library path (for extra shared libraries) */
    MCC_MagicExample_libpath_data,
    /* Number of directories in the load library path */
    0,

    /* MCR instance-specific runtime options */
    MCC_MagicExample_app_opts_data,
    /* Number of MCR instance-specific runtime options */
    0,

    /* MCR global runtime options */
    MCC_MagicExample_run_opts_data,
    /* Number of MCR global runtime options */
    0,
    
    /* Component preferences directory */
    "MagicExample_71C08FFE0B56F699EEDBC729EA4FA9E6",

    /* MCR warning status data */
    MCC_MagicExample_warning_state_data,
    /* Number of MCR warning status modifiers */
    1,

    /* Path to component - evaluated at runtime */
    NULL

};

#ifdef __cplusplus
}
#endif


