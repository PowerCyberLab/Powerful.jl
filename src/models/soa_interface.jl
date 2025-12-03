export to_struct_array

"""
Extract numerical fields from a vector model into a StructArray
"""
function to_struct_array(vec_data::T) where T
    # Get the base type name by removing "Vec" suffix
    type_str = string(nameof(typeof(vec_data)))
    @assert endswith(type_str, "Vec") "Input type must end with 'Vec'"
    
    base_type = Symbol(replace(type_str, "Vec" => ""))
    numerical_type = Symbol(string(base_type, "Numerical"))
    
    # Filter numerical fields using fieldnames and fieldtypes
    num_fields = filter(fn -> fieldtype(typeof(vec_data), fn).parameters[1] <: Number, 
                       fieldnames(typeof(vec_data)))
    
    # Create a NamedTuple of numerical fields
    nums = NamedTuple(fn => getfield(vec_data, fn) for fn in num_fields)
    
    # Create the StructArray with the constructed type
    return StructArray{eval(numerical_type)}(nums)
end

