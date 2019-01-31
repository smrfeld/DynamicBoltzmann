#include <string>
#include <vector>
#include <memory>
#include <map>

/************************************
* Namespace for bmla
************************************/

namespace bmla {
    
	/****************************************
	Filename
	****************************************/

	struct FName {

		std::string name;
		bool binary;

		FName(std::string name, bool binary);
	};

	/****************************************
	Filename collection
	****************************************/

	class FNameColl {

	private:

		// Filenames
		std::vector<FName> _fnames;
		std::vector<int> _idxs;

	public:

		/********************
		Get fnames
		********************/

		const std::vector<FName>& get_fnames_all() const;
		const FName& get_fname(int idx) const;
        std::vector<FName> get_fnames_batch(std::vector<int> idxs) const;

		/********************
		Add fname
		********************/

		void add_fname(FName fname);

		/********************
		Get random subset
		********************/

		std::vector<int> get_random_subset(int size);
	};
};
