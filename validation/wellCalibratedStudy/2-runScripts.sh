lphy_dir="$(cd ./lphyScripts && pwd)"
export JAVA_HOME=~/zulu25.32.21-ca-fx-jdk25.0.2-linux_x64
export PATH=$JAVA_HOME/bin:$PATH

for file in "${lphy_dir}"/*.lphy;
do
    cd ../../
    ~/apache-maven-3.9.14/bin/mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="convert -r 200 -l 30000000 ${file}"
#    use this when you want MRCAPrior instead of CalibrationPrior
#    ~/apache-maven-3.9.14/bin/mvn -pl calibratedcpp-lphybeast-launcher exec:exec -Dlphybeast.args="convert -r 200 -MRCAPrior -l 30000000 ${file}"
    cd - > /dev/null
done