#include <gifgen.h>

static char gdFontTinyData[] = {

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,1,1,1,0,
1,1,1,1,1,
0,1,1,1,0,
0,0,1,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,1,
0,1,0,1,0,
0,0,1,0,1,
0,1,0,1,0,
0,0,1,0,1,
0,1,0,1,0,
0,0,1,0,1,

0,1,0,1,0,
0,1,0,1,0,
0,1,1,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,1,1,1,
0,0,0,1,0,
0,0,0,1,0,

1,1,1,0,0,
1,0,0,0,0,
1,1,0,0,0,
1,0,1,1,1,
1,0,1,0,0,
0,0,1,1,0,
0,0,1,0,0,
0,0,1,0,0,

0,1,1,0,0,
1,0,0,0,0,
0,1,1,0,0,
0,0,0,0,0,
0,0,1,1,0,
0,0,1,0,1,
0,0,1,1,0,
0,0,1,0,1,

1,0,0,0,0,
1,0,0,0,0,
1,0,0,0,0,
1,1,1,0,0,
0,0,1,1,1,
0,0,1,0,0,
0,0,1,1,0,
0,0,1,0,0,

0,0,0,0,0,
0,1,1,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,1,1,1,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,1,1,1,0,
0,0,1,0,0,
0,0,0,0,0,
0,1,1,1,0,
0,0,0,0,0,
0,0,0,0,0,

1,0,0,1,0,
1,1,0,1,0,
1,0,1,1,0,
1,0,0,1,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,1,1,

1,0,1,0,0,
1,0,1,0,0,
1,0,1,0,0,
0,1,0,0,0,
0,0,1,1,1,
0,0,0,1,0,
0,0,0,1,0,
0,0,0,1,0,

0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
1,1,1,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
1,1,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,1,1,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,

0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,1,1,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
1,1,1,1,1,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,

0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,0,
0,0,0,0,0,

0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,1,1,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,

0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
1,1,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,

0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
1,1,1,1,1,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
1,1,1,1,1,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,

0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,

0,0,0,0,0,
0,0,0,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,0,1,0,0,
0,0,0,1,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,0,1,0,0,
0,0,0,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
1,1,1,1,1,
0,1,0,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,1,0,
1,1,1,1,1,
0,0,1,0,0,
1,1,1,1,1,
0,1,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,1,1,0,
0,1,0,0,1,
1,1,1,0,0,
0,1,0,0,0,
0,1,1,1,0,
1,1,0,1,1,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,0,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,1,0,1,0,
0,1,0,1,0,
1,1,1,1,1,
0,1,0,1,0,
1,1,1,1,1,
0,1,0,1,0,
0,1,0,1,0,
0,0,0,0,0,

0,0,1,0,0,
0,1,1,1,0,
1,0,1,0,0,
0,1,1,1,0,
0,0,1,0,1,
0,1,1,1,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,1,0,1,0,
0,0,1,0,0,
0,1,0,1,0,
0,0,0,1,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,1,0,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,1,0,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,1,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,1,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,1,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,1,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,0,1,1,0,
0,1,1,1,1,
0,0,1,1,0,
0,1,0,0,1,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,0,1,0,0,
1,1,1,1,1,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,0,0,
0,1,1,1,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,1,
0,0,0,0,1,
0,0,0,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,1,0,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,1,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,0,0,0,1,
0,0,1,1,0,
0,1,0,0,0,
0,1,1,1,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,0,0,1,0,
0,0,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,1,0,
0,0,1,1,0,
0,1,0,1,0,
0,1,1,1,1,
0,0,0,1,0,
0,0,0,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,1,
0,1,0,0,0,
0,1,1,1,0,
0,0,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,0,
0,1,0,1,0,
0,1,1,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,1,
0,0,0,1,0,
0,0,0,1,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,1,1,
0,0,1,0,1,
0,0,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,0,1,1,0,
0,0,0,0,0,
0,0,1,1,0,
0,0,1,1,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,0,1,1,0,
0,0,0,0,0,
0,0,1,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,0,1,0,0,
0,0,0,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,0,
0,0,0,0,0,
0,1,1,1,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,0,1,0,0,
0,0,0,1,0,
0,0,0,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,1,0,1,0,
0,0,0,1,0,
0,0,1,0,0,
0,0,0,0,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,1,1,0,
0,1,0,0,1,
1,0,0,1,1,
1,0,1,0,1,
1,0,1,0,1,
1,0,0,1,0,
0,1,0,0,0,
0,0,1,1,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,1,
0,1,0,0,0,
0,1,1,1,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,1,1,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,1,
0,1,0,0,0,
0,1,1,1,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,0,
0,1,0,1,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,1,
0,0,0,0,1,
0,0,0,0,1,
0,0,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,0,1,0,
0,1,1,0,0,
0,1,0,1,0,
0,1,0,1,0,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,1,1,1,
0,1,1,1,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,1,0,1,
0,1,1,1,1,
0,1,0,1,1,
0,1,0,1,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,0,
0,1,0,0,0,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,0,1,
0,1,0,1,1,
0,0,1,1,0,
0,0,0,0,1,

0,0,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,0,
0,1,0,1,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,0,1,0,0,
0,0,0,1,0,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
1,1,1,1,1,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,1,
0,1,1,1,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
1,0,0,0,1,
1,0,0,0,1,
0,1,0,1,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,1,
0,0,0,0,1,
0,0,0,1,0,
0,0,1,0,0,
0,1,0,0,0,
0,1,1,1,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,0,1,0,0,
0,0,0,1,0,
0,0,0,0,1,
0,0,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,1,0,
0,0,0,1,0,
0,0,0,1,0,
0,0,0,1,0,
0,0,0,1,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
1,1,1,1,1,

0,0,0,0,0,
0,1,1,0,0,
0,1,0,0,0,
0,0,1,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,0,1,
0,1,0,1,1,
0,1,0,1,1,
0,0,1,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,0,
0,1,0,0,0,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,1,
0,0,0,0,1,
0,0,1,0,1,
0,1,0,1,1,
0,1,0,1,1,
0,0,1,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,1,0,
0,1,1,1,1,
0,1,0,0,0,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,1,0,
0,0,1,0,1,
0,0,1,0,0,
0,1,1,1,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,0,1,1,1,
0,0,0,0,1,
0,0,1,1,0,

0,0,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,0,0,0,0,
0,1,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,1,0,
0,0,0,0,0,
0,0,0,1,0,
0,0,0,1,0,
0,0,0,1,0,
0,1,0,1,0,
0,0,1,0,0,

0,0,0,0,0,
0,1,0,0,0,
0,1,0,0,0,
0,1,0,0,1,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,1,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,0,1,0,
1,0,1,0,1,
1,0,1,0,1,
1,0,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,1,0,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,0,
0,1,0,0,1,
0,1,1,1,0,
0,1,0,0,0,
0,1,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,1,1,1,
0,1,0,0,1,
0,0,1,1,1,
0,0,0,0,1,
0,0,0,0,1,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,0,1,0,
0,1,1,0,1,
0,1,0,0,0,
0,1,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,0,
0,1,1,0,0,
0,0,0,1,0,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,1,1,1,0,
0,0,1,0,0,
0,0,1,0,1,
0,0,0,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,0,0,1,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,0,1,0,
0,1,0,1,0,
0,1,0,1,0,
0,0,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
1,0,0,0,1,
1,0,1,0,1,
1,0,1,0,1,
0,1,1,1,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,0,0,1,
0,0,1,1,0,
0,0,1,1,0,
0,1,0,0,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,0,0,1,
0,1,0,0,1,
0,0,1,1,1,
0,1,0,0,1,
0,0,1,1,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,1,1,1,1,
0,0,0,1,0,
0,0,1,0,0,
0,1,1,1,1,
0,0,0,0,0,

0,0,0,1,1,
0,0,1,0,0,
0,0,0,1,0,
0,1,1,0,0,
0,0,0,1,0,
0,0,1,0,0,
0,0,0,1,1,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,1,0,0,
0,0,0,0,0,

0,1,1,0,0,
0,0,0,1,0,
0,0,1,0,0,
0,0,0,1,1,
0,0,1,0,0,
0,0,0,1,0,
0,1,1,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,1,0,1,
0,1,0,1,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,

0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0,
0,0,0,0,0
};

static gdFont gdFontTinyRep = {
	128,
	0,
	5,
	8,
   1,
	gdFontTinyData
};

gdFontPtr gdFont5X8 = &gdFontTinyRep;
