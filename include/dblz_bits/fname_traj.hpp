#include <string>
#include <vector>
#include <memory>
#include <map>

/************************************
* Namespace for bmla
************************************/

namespace dblz {
    
	/****************************************
	Filename
	****************************************/

	struct FName {

		std::string name;
		bool binary;

		FName(std::string name, bool binary);
	};

    /****************************************
     Filename traj
     ****************************************/

    typedef std::vector<FName> FNameTraj;
    
	/****************************************
	Filename collection
	****************************************/

	class FNameTrajColl {

	private:

		// Filenames
		std::vector<FNameTraj> _fnames;
		std::vector<int> _idxs;

	public:

		/********************
		Get fnames
		********************/

		const std::vector<FNameTraj>& get_fname_trajs_all() const;
		const FNameTraj& get_fname_traj(int idx) const;
        std::vector<FNameTraj> get_fname_trajs_batch(std::vector<int> idxs) const;

		/********************
		Add fname
		********************/

		void add_fname_traj(FNameTraj fname_traj);

		/********************
		Get random subset
		********************/

		std::vector<int> get_random_subset(int size);
	};
};
