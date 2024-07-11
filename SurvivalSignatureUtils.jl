module SurvivalSignatureUtils

# ==============================================================================

# Function to ensure input array is treated as a column array
function ensure_column_array!(arr::Array)
    if size(arr, 1) > 1 && size(arr, 2) > 1
        transpose(arr)
    end
    return arr
end

# Function to ensure input array is treated as a row array
function ensure_row_array!(arr::Array)
    if size(arr, 1) > 1 && size(arr, 2) > 1
        arr[:] = transpose(arr)
    end
    return arr
end

function print_array(array::Array)
    array = ensure_row_array!(array)
    for row in eachrow(array)
        println(row)
    end
end

function header()
    open("Modules/header.txt", "r") do file
        # Read the entire contents of the file
        contents = read(file, String)
        # Print the contents to the console
        println(contents)
    end
end

# ================================ MACROS ======================================

export @examine

macro examine(var, value_bool::Bool=true)
    esc_var = esc(var)  # Escape the variable `var` to ensure it's captured correctly
    var_name = string(var)  # Get the variable name as a string

    quote
        println("───────────────────────────────────────")
        println("Variable name: ", $var_name)  # Print the variable name directly
        println("───────────────────────────────────────")

        if $value_bool
            println("Value: ", $esc_var)  # Print the value directly
        end

        if isa($esc_var, AbstractArray)
            println("Size: ", size($esc_var))  # Calculate size directly
        end

        println("Type: ", typeof($esc_var))  # Get type directly

        println("───────────────────────────────────────")
    end
end

export @allow_nothing
macro allow_nothing(expr)

    # THE EXISTENCE OF COMMENTS WITHIN THE struct BREAKS THIS MACRO #

    # Ensure the input is a struct definition
    if expr.head != :struct
        error("@allow_nothing must be used with a struct definition")
    end

    # Extract the struct name and fields block
    struct_name = expr.args[2]
    fields_block = expr.args[3]

    # Ensure fields_block is a quoted expression
    if fields_block.head != :block
        error("Unexpected struct definition format")
    end

    # Process each field to include Union{Nothing, ...} or set to Any if no type annotation
    new_fields = map(
        field -> begin
            if field isa Expr && field.head == :(::)  # Check if the field has a type annotation
                field_name = field.args[1]
                field_type = field.args[2]
                :($field_name::Union{Nothing,$field_type})
            elseif field isa LineNumberNode  # Ignore line number nodes (comments)
                nothing
            else  # If no type annotation, set the field type to Any
                :($field::Union{Nothing,Any})
            end
        end, fields_block.args
    )

    # Construct the new struct definition
    return esc(
        quote
            mutable struct $struct_name
                $(new_fields...)
            end
        end,
    )
end
# ==============================================================================

export *

end
