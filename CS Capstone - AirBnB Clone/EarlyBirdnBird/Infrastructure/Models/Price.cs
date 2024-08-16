using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.ComponentModel.DataAnnotations.Schema;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Infrastructure.Models
{
    public class Price
    {
        [Key]
        public int Id { get; set; }

        [ForeignKey("RentalId")]
        public int RentalId { get; set; }


        [DataType(DataType.Date)]
        public DateTime PriceDateStart { get; set; }

        [DataType(DataType.Date)]
        public DateTime PriceDateEnd { get; set; }

        public double? PricePerNight { get; set; }
    }
}
