module XYZReader

export read_xyz_ascii

"""
    read_xyz_ascii(path; comment = '#')

Read a simple ASCII XYZ file and return coordinate vectors `(x, y, z)`.

Accepted formats per line (whitespace or comma separated):
- `x y z`
- `x, y, z`
- `x y z <extra...>`  (extra columns are ignored)

Blank lines and lines starting with `comment` are skipped.
All values are parsed as `Float64`.
"""
function read_xyz_ascii(path::AbstractString; comment::AbstractString = "#")
    xs = Float64[]; ys = Float64[]; zs = Float64[]
    open(path, "r") do io
        while !eof(io)
            line = strip(readline(io))
            isempty(line) && continue
            startswith(line, comment) && continue
            # allow commas or whitespace
            line = replace(line, ',' => ' ')
            toks = split(line)
            length(toks) < 3 && continue
            try
                x = parse(Float64, toks[1])
                y = parse(Float64, toks[2])
                z = parse(Float64, toks[3])
                push!(xs, x); push!(ys, y); push!(zs, z)
            catch
                # Ignore malformed lines
            end
        end
    end
    return xs, ys, zs
end

end # module
