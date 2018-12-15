//
// Created by morfeush22 on 15.12.18.
//

#include "base.h"

void basicTest() {
    ASSERT_EQUAL(1, 1);
    struct Data data = loadDataFromFile("sources/simple_eq");
    printData(data);
}

int main(int argc, char** argv) {
    basicTest();
}
