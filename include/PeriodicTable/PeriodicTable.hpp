#pragma once

#include <string>
#include <unordered_map>

struct Element {
    std::string name;
    int z;
    double mass;
    int n_e;
    int n_v_e;
};


class PeriodicTable {

    private:
        std::unordered_map<std::string, Element> elements_by_symbol;
        std::unordered_map<int, Element> elements_by_atomic_number;

    public:

        PeriodicTable() {

             Element HYDROGEN = {"H", 1, 1.008, 1, 1};
             Element CARBON = {"C", 6, 12.011, 6, 4};
             Element NITROGEN = {"N", 7, 14.007, 7, 5};
             Element OXYGEN = {"O", 8, 15.999, 8, 6};
             Element FLUORINE = {"F", 9, 18.998, 9, 7};

            add_element(HYDROGEN);
            add_element(CARBON);
            add_element(NITROGEN);
            add_element(OXYGEN);
            add_element(FLUORINE);
        }

        void add_element(Element& e) {
            elements_by_symbol[e.name] = e;
            elements_by_atomic_number[e.z] = e;
        }

        Element get_element(std::string name) {
            auto it = this->elements_by_symbol.find(name);
            if (it != elements_by_symbol.end()) {
                return it->second;
            } else {
                std::cout << "ELEMENT " << name << " NOT IN TABLE!!" << std::endl;
                exit(1);
            }
        }
        Element get_element(int z) {
            auto it = this->elements_by_atomic_number.find(z);
            if (it != elements_by_atomic_number.end()) {
                return it->second;
            } else {
                std::cout << "ELEMENT " << z << " NOT IN TABLE!!" << std::endl;
                exit(1);
            }
        }
};


