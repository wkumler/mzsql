using SQLite, DataFrames, StatsPlots

db = SQLite.DB("chrom_data/chrom_data.sqlite")
targets = [147.0764, 148.0734, 149.0705]
window(m) = (m*(1 - 2/1e6), m*(1 + 2/1e6))
clause = join(["(mz BETWEEN $(mzmin) AND $(mzmax))" for (mzmin, mzmax) in window.(targets)], " OR ")
query = "SELECT * FROM MS1 WHERE ($clause) AND (rt BETWEEN 10 AND 12)"
df = DBInterface.execute(db, query) |> DataFrame
SQLite.close(db)

p1 = @df df[round.(df.mz) .== 147, :] plot(:rt, :int, group=:filename, color="black", title="Unlabeled")
p2 = @df df[round.(df.mz) .== 148, :] plot(:rt, :int, group=:filename, color="#0b505c", title="15Nx1")
p3 = @df df[round.(df.mz) .== 149, :] plot(:rt, :int, group=:filename, color="#028e34", title="15Nx2")
plot(p1, p2, p3, layout = (3,1), dpi=600)


using LaTeXStrings, Plots
default(fontfamily="Computer Modern")
p1 = @df df[round.(df.mz) .== 147, :] plot(:rt, :int, group=:filename, color="black", legend=false, 
title="Unlabeled", xlims=(10, 12), titlefontsize = 10)
p2 = @df df[round.(df.mz) .== 148, :] plot(:rt, :int, group=:filename, color="#0b505c", legend=false, 
title=L"\textrm{^{15}N_1\ isotope}", xlims=(10, 12), titlefontsize = 10)
p3 = @df df[round.(df.mz) .== 149, :] plot(:rt, :int, group=:filename, color="#028e34", legend=false, 
title=L"\textrm{^{15}N_2\ isotope}", xlims=(10, 12), titlefontsize = 10)
plot(p1, p2, p3, layout = (3,1), dpi=600, size=(600, 400))
savefig("chrom_sqlite_julia_fig.png")
