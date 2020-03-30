#include <iostream>

#include "OpDyn.hh"

using namespace Utopia::Models::OpDyn;


int main (int argc, char** argv)
{
    try {
        // Initialize the PseudoParent from config file path
        Utopia::PseudoParent pp(argv[1]);

        auto model_cfg = pp.get_cfg()["OpDyn"];
        auto ageing = Utopia::get_as<std::string>("user_ageing", model_cfg);
        auto media = Utopia::get_as<std::string>("media_status", model_cfg);
        
        if (ageing=="on") {
            if (media=="on") {
                OpDyn<Ageing_and_Media> model("OpDyn", pp);
                model.run();
            }
            else if (media=="off") {
                OpDyn<Ageing> model("OpDyn", pp);
                model.run();

            }
            else {
                throw std::invalid_argument("Media mode {} unknown! Set media "
                                            "to either 'on' or 'off'");
            }
        }
        
        else if (ageing=="off") {
            if (media=="on") {
                OpDyn<Media> model ("OpDyn", pp);
                model.run();

            }
            if (media=="off") {
                OpDyn<None> model("OpDyn", pp);
                model.run();

            }
            else {
                throw std::invalid_argument("Media mode {} unknown! Set media "
                                            "to either 'on' or 'off'");
            }
        }
        
        else {
            throw std::invalid_argument("Ageing mode {} unkown! Set ageing "
                                        " to either 'on' or 'off'");
        }
        return 0;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Exception occured!" << std::endl;
        return 1;
    }
}
