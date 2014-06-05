#ifndef GRID_H
#define GRID_H
float GridGetValue(grid *G, float x, float y, float z);
float GridGetValueIndex(grid *G, int index);
grid * GridCreate(void);
void GridNullify(grid *C);
void GridFree(grid *C);
#endif
