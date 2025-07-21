#pragma once
#include <vector>
#include <memory>
#include <string>
#include "BaseModel.h"

class Models {
public:
    void registerModel(std::shared_ptr<BaseModel> model);
    std::shared_ptr<BaseModel> getModel(const std::string& name) const;
    std::vector<std::string> listModelNames() const;

private:
    std::vector<std::shared_ptr<BaseModel>> modelList;
};
