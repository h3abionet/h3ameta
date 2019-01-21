

params.message = "hello"

data_ch = Channel.from(params.message)

process banner {
  input:
    val(message) from data_ch
  output:
    file(result) into banner_ch
  script:
    result = "banner.txt"
    """
    banner ${params.message} >  $result
    """
}


process countRows {
  validExitStatus 0, 1
  input:
    file(result) from banner_ch
  each col from 1..80
  output:
    stdout col_ch
  script:
  """
    cut -b $col $result > col.txt
    grep -c \\#  col.txt
  """
}

process combineResults {
  input:
  val(res) from col_ch.toList()
  output:
    set  val(size), file(result) into comb_ch
    val(res) into another_ch
  script:
    size = res.size()
    data = res*.trim().join(" ")
    """
     count_num.py $data > result
   """
}

inp_file_ch = Channel.fromPath(params.compfile)

process compareCh {
  input:
    set val(size), file(result) from comb_ch
    file(comp) from inp_file_ch
  publishDir 'result'
  output:
    file("output.txt")
  script:
  """
    num_full=`cat $result`
    num_comp=`wc -l $comp | cut -f 1 -d\\ `
    if [ \$num_full -gt \$num_comp ]; then
       echo "Columns" > output.txt
    else
       echo "File $result" > output.txt
    fi
  """
}
 
