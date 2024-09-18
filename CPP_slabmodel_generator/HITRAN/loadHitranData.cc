#include "Hitran.ih"

void Hitran::loadHitranData(double lowerLam, double upperLam) {
    readFixedWidthFile();
    filterData(lowerLam, upperLam);
}