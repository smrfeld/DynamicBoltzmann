/************************************
* Namespace for bmla
************************************/

namespace bmla {

    enum class Solver : unsigned int { SGD, NESTEROV, ADAM };

    enum class CDModeAsleep : unsigned int { PERSISTENT_CD, START_FROM_DATA, START_FROM_RANDOM };
        
    struct OptionsAsleepPersistentCD {
                
        OptionsAsleepPersistentCD();
        ~OptionsAsleepPersistentCD();
    };
    
    struct OptionsAsleepStartFromRandom {
        
        // Init random - binary?
        bool start_from_binary_visible = true;
        //    Binary hidden?
        bool start_from_binary_hidden = true;
        
        OptionsAsleepStartFromRandom();
        ~OptionsAsleepStartFromRandom();
    };
    
    struct OptionsAsleepStartFromData {
        // Start from binary?
        bool start_from_binary_hidden = true;
        
        OptionsAsleepStartFromData();
        ~OptionsAsleepStartFromData();
    };

};
