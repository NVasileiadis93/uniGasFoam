template<class Type>
const Foam::fvMesh& Foam::CloudWithModels<Type>::mesh() const
{
    return mesh_;
}


template<class Type>
const Foam::IOdictionary&
Foam::CloudWithModels<Type>::particleProperties() const
{
    return particleProperties_;
}


template<class Type>
const Foam::IOdictionary& Foam::CloudWithModels<Type>::outputProperties() const
{
    return outputProperties_;
}


template<class Type>
Foam::IOdictionary& Foam::CloudWithModels<Type>::outputProperties()
{
    return outputProperties_;
}


template<class Type>
const Foam::cloudSolution& Foam::CloudWithModels<Type>::solution() const
{
    return solution_;
}


template<class Type>
inline Foam::List<Foam::DynamicList<Type*>>&
Foam::CloudWithModels<Type>::cellOccupancy()
{
    if (!cellOccupancyPtr_)
    {
        buildCellOccupancy();
    }

    return *cellOccupancyPtr_;
}
